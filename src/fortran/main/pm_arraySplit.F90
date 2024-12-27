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
!>  This module contains procedures and generic interfaces for splitting arrays of various types at the specified instances of occurrence of pattern.
!>
!>  \benchmarks
!>
!>  \benchmark{scalarSep_vs_vectorSep, The runtime performance of [setSplit](@ref pm_arraySplit::setSplit) for scalar vs. vector input `sep` argument.}
!>  \include{lineno} benchmark/pm_arraySplit/scalarSep_vs_vectorSep/main.F90
!>  \compilefb{scalarSep_vs_vectorSep}
!>  \postprocb{scalarSep_vs_vectorSep}
!>  \include{lineno} benchmark/pm_arraySplit/scalarSep_vs_vectorSep/main.py
!>  \visb{scalarSep_vs_vectorSep}
!>  \image html benchmark/pm_arraySplit/scalarSep_vs_vectorSep/benchmark.scalarSep_vs_vectorSep.runtime.png width=1000
!>  \image html benchmark/pm_arraySplit/scalarSep_vs_vectorSep/benchmark.scalarSep_vs_vectorSep.runtime.ratio.png width=1000
!>  \moralb{scalarSep_vs_vectorSep}
!>      -#  The procedures under the generic interface [setSplit](@ref pm_arraySplit::setSplit) take both scalar and vector `sep` arguments.<br>
!>          As evidenced by the above benchmark, when the input `sep` is vector of length `1`, it is much faster, **up to 4X**,
!>          to pass `sep` as a scalar instead of a whole array of length `1`.<br>
!>          Note that this benchmark is likely irrelevant to removing substrings from Fortran strings.<br>
!>
!>  \test
!>  [test_pm_arraySplit](@ref test_pm_arraySplit)
!>
!>  \todo
!>  \pmed
!>  A benchmark comparing the performance of output index array vs. output jagged array would be informative here.<br>
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

module pm_arraySplit

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_arraySplit"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the parts of the input array split at the requested occurrences of the input `sep`.
    !>
    !>  \details
    !>  If an input vector of `instance` is specified, containing the indices of the specific instances of `sep` at
    !>  which the `array` must be split, then the input `array` will be split only at those specific requested instances.
    !>
    !>  \param[out] field       :   The output object that can be any of the following:<br>
    !>                              <ol>
    !>                                  <li>    The output `allocatable` array of shape `(2, :)` of type `integer` of default kind \IK, <br>
    !>                                          whose length along the second dimension is the number of splits of the input array.<br>
    !>                                          The two elements along the first dimension represent the starting and the ending indices of each split part in `array` such that,<br>
    !>                                          <ol>
    !>                                              <li>    The subset `field(1, :)` contains the indices of the beginnings of the split parts of the input `array`.<br>
    !>                                              <li>    The subset `field(2, :)` contains the indices of the end points of the split parts of the input `array`.<br>
    !>                                          </ol>
    !>                                  <li>    The output `allocatable` derived object (jagged-array) of shape `(:)` of PDT type of either, <br>
    !>                                          <ol>
    !>                                              <li>    [css_pdt](@ref pm_container::css_pdt) of the same kind as the input `array` argument, or<br>
    !>                                              <li>    [cvs_pdt](@ref pm_container::cvs_pdt) of the same kind as the input `array` argument, or<br>
    !>                                              <li>    [cvi_pdt](@ref pm_container::cvi_pdt) of the same kind as the input `array` argument, or<br>
    !>                                              <li>    [cvl_pdt](@ref pm_container::cvl_pdt) of the same kind as the input `array` argument, or<br>
    !>                                              <li>    [cvc_pdt](@ref pm_container::cvc_pdt) of the same kind as the input `array` argument, or<br>
    !>                                              <li>    [cvr_pdt](@ref pm_container::cvr_pdt) of the same kind as the input `array` argument<br>
    !>                                          </ol>
    !>                                  <li>    The output `allocatable` derived object (jagged-array) of shape `(:)` of type of either, <br>
    !>                                          <ol>
    !>                                              <li>    [css_type](@ref pm_container::css_type) of the same **default** kind \SK as the input `array` argument, or<br>
    !>                                              <li>    [cvs_type](@ref pm_container::cvs_type) of the same **default** kind \SK as the input `array` argument, or<br>
    !>                                              <li>    [cvi_type](@ref pm_container::cvi_type) of the same **default** kind \IK as the input `array` argument, or<br>
    !>                                              <li>    [cvl_type](@ref pm_container::cvl_type) of the same **default** kind \LK as the input `array` argument, or<br>
    !>                                              <li>    [cvc_type](@ref pm_container::cvc_type) of the same **default** kind \CK as the input `array` argument, or<br>
    !>                                              <li>    [cvr_type](@ref pm_container::cvr_type) of the same **default** kind \RK as the input `array` argument<br>
    !>                                          </ol>
    !>                                          containing pieces of the input `array` split at the requested instances of the input `sep`.<br>
    !>                                          Note that the vector contained within the container class (,that is, the allocatable component of the jagged-array)
    !>                                          must have the same type and kind as the input `array`.
    !>                              </ol>
    !>  \param[in]  array       :   The input `contiguous` array of shape `(:)` of the same type and kind as the `Value` component of the output argument `field`, which can be either <br>
    !>                              <ol>
    !>                                  <li>    type `character` of kind \SKALL, or <br>
    !>                                  <li>    type `logical` of kind \LKALL, or <br>
    !>                                  <li>    type `integer` of kind \IKALL, or <br>
    !>                                  <li>    type `complex` of kind \CKALL, or <br>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ol>
    !>                              or,
    !>                              <ol>
    !>                                  <li>    a scalar assumed-length `character` of kind \SKALL, <br>
    !>                              </ol>
    !>                              that will be split at the specified instances of the input `sep`.
    !>  \param[in]  sep         :   The input `contiguous` array of shape `(:)` or scalar of the same type and kind as the input `array`,
    !>                              containing the pattern at which the input `array` will be split.
    !>  \param      iseq        :   The `external` user-specified function that takes either two input assumed-length `character` arguments (if
    !>                              if the input `array` is also an assumed-length `character`) or two array-valued **explicit-shape** arguments of the same
    !>                              type and kind as the input `array`.<br>
    !>                              It must return a scalar of type `logical` of default kind \LK that is `.true.` if all elements of
    !>                              the two input arguments are equivalent (e.g., equal) according to the user-defined criterion, otherwise, it is `.false.`.<br>
    !>                              The the input `array` is array-valued, then the last argument to `iseq` is the length of the input `sep`.<br>
    !>                              This user-defined equivalence check is extremely useful where an equivalence test other than exact identity is needed,
    !>                              for example, when the array segments should match the input `sep` only within a given threshold or,
    !>                              when the case-sensitivity in character comparisons do not matter.<br>
    !>                              In such cases, user can define a custom equivalence criterion within the user-defined external function `iseq` to achieve the goal.<br>
    !>                              The following illustrates the generic interface of `iseq` when `sep` is array-valued,
    !>                              \code{.F90}
    !>                                  function iseq(Segment, sep, sepSize) result(equivalent)
    !>                                      use pm_kind, only: IK, LK
    !>                                      integer(IK)     , intent(in)    :: sepSize
    !>                                      TYPE(KIND)      , intent(in)    :: Segment(sepSize), sep(sepSize)
    !>                                      logical(LK)                     :: equivalent
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array` and `sepSize = size(sep, kind = IK)`.<br>
    !>                              The following illustrates the generic interface of `iseq` when `sep` is scalar-valued,
    !>                              \code{.F90}
    !>                                  function iseq(segment, sep) result(equivalent)
    !>                                      use pm_kind, only: IK, LK
    !>                                      TYPE(KIND)      , intent(in)    :: segment, sep
    !>                                      logical(LK)                     :: equivalent
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type (`character`) and kind of the input argument `array`.<br>
    !>                              (**optional**, the default equivalence operator is `.eqv.` if the input `array` is `logical`, otherwise `==`).
    !>  \param[in]  instance    :   The input `contiguous` array of shape `(:)` of type `integer` of default kind \IK,
    !>                              containing the instances of the input `sep` in the input `array` at which the `array` should be split.<br>
    !>                              Any element of `instance` that points to an out-of-scope instance of `sep` in the input `array` will be ignored.<br>
    !>                              Any element of `instance` that is negatively valued will be counted from end of the input `array`.<br>
    !>                              For example, `instance = [2,-1]` requests splitting at the second instance of `sep` in `array` from the beginning and
    !>                              splitting at the first instance of `sep` starting from the end of `array`.<br>
    !>                              (**optional**, the default value corresponds to splitting at all instances of `sep` in `array`)
    !>  \param[in]  sorted      :   The input scalar `logical` of default kind \LK indicating whether the elements of the specified input `instance` are all in ascending-order.<br>
    !>                              This includes the negative elements of `instance` **after** they are translated to the corresponding **positive** instances from
    !>                              the beginning of the input `array`.<br>Setting `sorted = .true.` will lead to faster runtime of the procedure.<br>
    !>                              However, the onus will be strictly on the user to ensure all elements of `instance` are in ascending-order.<br>
    !>                              This is generally not an easy guarantee to make if there are negative elements in `instance`.<br>
    !>                              Therefore, set `sorted = .true.` <b>if and only if</b> you can guarantee the validity of the condition.<br>
    !>                              (**optional**, default = `.false.`. It can be present as input argument <b>if and only if</b> the input argument `instance` is present.)
    !>  \param[in]  unique      :   The input scalar `logical` of default kind \LK indicating whether the elements of the specified input `instance` are all unique.<br>
    !>                              This includes the negative elements of `instance` **after** they are translated to the corresponding **positive** instances from
    !>                              the beginning of the input `array`.<br>Setting `unique = .true.` will lead to faster runtime of the procedure.<br>
    !>                              However, the onus will be strictly on the user to ensure all elements of `instance` are unique.<br>
    !>                              This is generally not an easy guarantee to make if there are negative elements in `instance`.<br>
    !>                              Therefore, set `unique = .true.` <b>if and only if</b> you can guarantee the validity of the condition.<br>
    !>                              (**optional**, default = `.false.`. It can be present as input argument <b>if and only if</b> the input argument `instance` is present.)
    !>  \param[in]  keep        :   The input scalar `logical` of default kind \LK. If `.true.`, all instances of `sep` will be kept in the output as if each one results from splitting.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_arraySplit, only: setSplit
    !>
    !>      call setSplit(field, array, sep, keep = keep)
    !>      call setSplit(field, array, sep, iseq, keep = keep)
    !>      call setSplit(field, array, sep, instance, sorted = sorted, unique = unique, keep = keep)
    !>      call setSplit(field, array, sep, iseq, instance, sorted = sorted, unique = unique, keep = keep)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \warning
    !>  The procedures under this generic interface are `impure` when the user-specified `external` procedure `iseq` is specified as input argument.<br>
    !>
    !>  \warning
    !>  Note that in Fortran, trailing blanks are ignored in character comparison, that is, `"Fortran" == "Fortran "` yields `.true.`.<br>
    !>
    !>  \see
    !>  [setReplaced](@ref pm_arrayReplace::setReplaced)<br>
    !>  [getReplaced](@ref pm_arrayReplace::getReplaced)<br>
    !>  [getInserted](@ref pm_arrayInsert::getInserted)<br>
    !>  [setInserted](@ref pm_arrayInsert::setInserted)<br>
    !>  [getRemoved](@ref pm_arrayRemove::getRemoved)<br>
    !>  [setRemoved](@ref pm_arrayRemove::setRemoved)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_arraySplit/setSplit/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_arraySplit/setSplit/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arraySplit](@ref test_pm_arraySplit)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifort{2021.5}
    !>  \desc
    !>  The \ifort{2021.5} on WSL platform cannot correctly
    !>  set the PDT type alias in a module `use` statement.<br>
    !>  The following is a sample code demonstrating the issue,
    !>  \code{.F90}
    !>
    !>      module container_module
    !>
    !>          use iso_fortran_env
    !>
    !>          integer, parameter :: lkc = logical_kinds(1)
    !>          integer, parameter :: ikc = integer_kinds(1)
    !>
    !>          type :: container_ik(kind)
    !>              integer         , kind  :: kind = integer_kinds(1)
    !>              integer(kind)           :: value
    !>          end type
    !>
    !>          type :: container_lk(kind)
    !>              integer         , kind  :: kind = logical_kinds(1)
    !>              logical(kind)           :: value
    !>          end type
    !>
    !>      end module
    !>
    !>      module disp_module
    !>
    !>          use container_module
    !>
    !>          type :: disp_type
    !>          end type
    !>
    !>          interface show
    !>          module subroutine show_logical(disp, object)
    !>              use container_module, only: container => container_lk
    !>              class(disp_type)        , intent(inout)               :: disp
    !>             !type(container_lk(lkc)) , intent(in)                  :: object
    !>              type(container(lkc))    , intent(in)                  :: object
    !>          end subroutine
    !>          module subroutine show_integer(disp, object)
    !>              use container_module, only: container => container_ik
    !>              class(disp_type)        , intent(inout)               :: disp
    !>             !type(container_ik(ikc)) , intent(in)                  :: object
    !>              type(container(ikc))    , intent(in)                  :: object
    !>          end subroutine
    !>          end interface
    !>
    !>      end module
    !>
    !>      end
    !>
    !>  \endcode
    !>  The ifort{2021.5} generates the following error with the above code,
    !>  \code{bash}
    !>      error #5286: Ambiguous generic interface SHOW: previously declared specific procedure SHOW_LOGICAL is not distinguishable from this declaration. [SHOW_INTEGER]
    !>         !module subroutine show_integer(disp, object)
    !>      ----------------------^
    !>  \endcode
    !>  The error however is incorrect because the two interfaces are different.<br>
    !>  The issue appears to be with the renaming of the module entities within the interfaces to the same aliases ("container").<br>
    !>  If we uncomment the commented lines to use the type names explicitly, then the error is resolved.<br>
    !>  \remedy
    !>  For now, all PDT name aliases have been removed to bypass the \ifort bug.<br>
    !>
    !>  \todo
    !>  \plow For now, all PDT name aliases have been removed to bypass the \ifort bug.<br>
    !>  This can be reverted back to the original aliasing approach in the future once the Intel ifort bug is resolved.<br>
    !>  However, the whole point of using aliases was to make the development of the code easier and nicer.<br>
    !>  With the use of explicit names, there might really be no point in reviving the aliases other than code aesthetics.<br>
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to 2D input objects.
    !>
    !>  \todo
    !>  \pmed A `backward` optional argument can be added in the future to search for the input `sep` in the backward direction.
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    !                                   <li>    An output `contiguous` vector of shape `(:)` of type `integer` of default kind \IK, <br>
    !                                           representing the positions of the last elements of the split parts in the input `array` such that,<br>
    !                                           <ol>
    !                                               <li>    The subset `array(1 : field(1) - lenDelim)` represents the first
    !                                               <li>    The subset `field(2, :)` contains the indices of the end points of the split parts of the input `array`.<br>
    !                                           </ol>
    !                                           The two elements along the first dimension represent the starting and the ending indices of each split part in `array`.<br>
    !                                           The actual final number of splits of the input array will be output in the output argument `nsplit` such that,<br>
    !                                           <ol>
    !                                               <li>    The subset `field(1, 1 : nsplit)` contains the indices of the beginnings of the split parts of the input `array`.<br>
    !                                               <li>    The subset `field(2, 1 : nsplit)` contains the indices of the end points of the split parts of the input `array`.<br>
    !                                           </ol>
    !                                           This `contiguous` form of `field` is available **if and only if** the output argument `nsplit` is also present.<br>
    !                                           This `contiguous` form of `field` is particularly useful and performant by avoiding repeated allocations of the output `field`
    !                                           when many records with roughly equal number of split parts are to be split using this generic interface.<br>
    !                                           The only caveat is that the user must set the size of the second dimension of `field`
    !                                           to be larger than or equal to the actual number of splits that the algorithm may find.<br>
    !                                           This limitation can be readily overcome by pre-allocating `field` to theoretical maximum possible
    !                                           value for the output `nsplit`, which yields the shape `(1 : 2, 1 : 1 + lenArray / min(1, lenDelim))`,
    !                                           where `lenArray` is the length/size of the input `array` and `lenDelim` is length/size of the input `sep`.<br>
    !   \param[in]  nsplit      :   The output scalar of type `integer` of default kind \IK containing the number of splits identified in the input `array`,
    !                               such that `field(1 : nsplit)` represents the vector of positions of the last elements of split parts in the input `array`.<br>
    !                               (**optional**. It must be present *if and only** the output `field` is a pre-sized (allocated) **vector** of type `integer`.)

!    ! fixed-size index
!
!    interface setSplit
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D0_D0_SK5(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D0_D0_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)                    :: array
!        character(*,SKG)            , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D0_D0_SK4(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D0_D0_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)                    :: array
!        character(*,SKG)            , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D0_D0_SK3(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D0_D0_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)                    :: array
!        character(*,SKG)            , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D0_D0_SK2(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D0_D0_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)                    :: array
!        character(*,SKG)            , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D0_D0_SK1(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D0_D0_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)                    :: array
!        character(*,SKG)            , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D0_SK5(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D0_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D0_SK4(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D0_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D0_SK3(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D0_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D0_SK2(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D0_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D0_SK1(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D0_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D0_IK5(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D0_IK5
!#endif
!        use pm_kind, only: IKG => IK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D0_IK4(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D0_IK4
!#endif
!        use pm_kind, only: IKG => IK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D0_IK3(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D0_IK3
!#endif
!        use pm_kind, only: IKG => IK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D0_IK2(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D0_IK2
!#endif
!        use pm_kind, only: IKG => IK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D0_IK1(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D0_IK1
!#endif
!        use pm_kind, only: IKG => IK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if LK5_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D0_LK5(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D0_LK5
!#endif
!        use pm_kind, only: LKG => LK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK4_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D0_LK4(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D0_LK4
!#endif
!        use pm_kind, only: LKG => LK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK3_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D0_LK3(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D0_LK3
!#endif
!        use pm_kind, only: LKG => LK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK2_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D0_LK2(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D0_LK2
!#endif
!        use pm_kind, only: LKG => LK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK1_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D0_LK1(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D0_LK1
!#endif
!        use pm_kind, only: LKG => LK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if CK5_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D0_CK5(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D0_CK5
!#endif
!        use pm_kind, only: CKG => CK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK4_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D0_CK4(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D0_CK4
!#endif
!        use pm_kind, only: CKG => CK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK3_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D0_CK3(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D0_CK3
!#endif
!        use pm_kind, only: CKG => CK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK2_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D0_CK2(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D0_CK2
!#endif
!        use pm_kind, only: CKG => CK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK1_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D0_CK1(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D0_CK1
!#endif
!        use pm_kind, only: CKG => CK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D0_RK5(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D0_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D0_RK4(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D0_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D0_RK3(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D0_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D0_RK2(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D0_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D0_RK1(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D0_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)                    :: sep
!        logical(LK)                 , intent(in)    , optional      :: keep
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
!    module subroutine setSplitFixCusComDefIns_D0_D0_SK5(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D0_D0_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)                    :: array
!        character(*,SKG)            , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    module subroutine setSplitFixCusComDefIns_D0_D0_SK4(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D0_D0_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)                    :: array
!        character(*,SKG)            , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    module subroutine setSplitFixCusComDefIns_D0_D0_SK3(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D0_D0_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)                    :: array
!        character(*,SKG)            , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    module subroutine setSplitFixCusComDefIns_D0_D0_SK2(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D0_D0_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)                    :: array
!        character(*,SKG)            , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    module subroutine setSplitFixCusComDefIns_D0_D0_SK1(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D0_D0_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)                    :: array
!        character(*,SKG)            , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D0_SK5(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D0_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D0_SK4(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D0_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D0_SK3(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D0_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D0_SK2(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D0_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D0_SK1(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D0_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D0_IK5(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D0_IK5
!#endif
!        use pm_kind, only: IKG => IK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D0_IK4(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D0_IK4
!#endif
!        use pm_kind, only: IKG => IK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D0_IK3(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D0_IK3
!#endif
!        use pm_kind, only: IKG => IK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D0_IK2(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D0_IK2
!#endif
!        use pm_kind, only: IKG => IK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D0_IK1(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D0_IK1
!#endif
!        use pm_kind, only: IKG => IK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if LK5_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D0_LK5(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D0_LK5
!#endif
!        use pm_kind, only: LKG => LK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK4_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D0_LK4(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D0_LK4
!#endif
!        use pm_kind, only: LKG => LK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK3_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D0_LK3(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D0_LK3
!#endif
!        use pm_kind, only: LKG => LK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK2_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D0_LK2(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D0_LK2
!#endif
!        use pm_kind, only: LKG => LK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK1_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D0_LK1(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D0_LK1
!#endif
!        use pm_kind, only: LKG => LK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if CK5_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D0_CK5(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D0_CK5
!#endif
!        use pm_kind, only: CKG => CK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK4_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D0_CK4(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D0_CK4
!#endif
!        use pm_kind, only: CKG => CK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK3_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D0_CK3(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D0_CK3
!#endif
!        use pm_kind, only: CKG => CK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK2_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D0_CK2(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D0_CK2
!#endif
!        use pm_kind, only: CKG => CK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK1_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D0_CK1(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D0_CK1
!#endif
!        use pm_kind, only: CKG => CK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D0_RK5(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D0_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D0_RK4(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D0_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D0_RK3(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D0_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D0_RK2(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D0_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D0_RK1(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D0_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
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
!    PURE module subroutine setSplitFixDefComCusIns_D0_D0_SK5(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D0_D0_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)                    :: array
!        character(*,SKG)            , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D0_D0_SK4(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D0_D0_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)                    :: array
!        character(*,SKG)            , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D0_D0_SK3(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D0_D0_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)                    :: array
!        character(*,SKG)            , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D0_D0_SK2(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D0_D0_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)                    :: array
!        character(*,SKG)            , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D0_D0_SK1(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D0_D0_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)                    :: array
!        character(*,SKG)            , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D0_SK5(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D0_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D0_SK4(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D0_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D0_SK3(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D0_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D0_SK2(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D0_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D0_SK1(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D0_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D0_IK5(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D0_IK5
!#endif
!        use pm_kind, only: IKG => IK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D0_IK4(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D0_IK4
!#endif
!        use pm_kind, only: IKG => IK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D0_IK3(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D0_IK3
!#endif
!        use pm_kind, only: IKG => IK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D0_IK2(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D0_IK2
!#endif
!        use pm_kind, only: IKG => IK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D0_IK1(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D0_IK1
!#endif
!        use pm_kind, only: IKG => IK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if LK5_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D0_LK5(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D0_LK5
!#endif
!        use pm_kind, only: LKG => LK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK4_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D0_LK4(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D0_LK4
!#endif
!        use pm_kind, only: LKG => LK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK3_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D0_LK3(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D0_LK3
!#endif
!        use pm_kind, only: LKG => LK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK2_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D0_LK2(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D0_LK2
!#endif
!        use pm_kind, only: LKG => LK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK1_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D0_LK1(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D0_LK1
!#endif
!        use pm_kind, only: LKG => LK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if CK5_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D0_CK5(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D0_CK5
!#endif
!        use pm_kind, only: CKG => CK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK4_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D0_CK4(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D0_CK4
!#endif
!        use pm_kind, only: CKG => CK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK3_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D0_CK3(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D0_CK3
!#endif
!        use pm_kind, only: CKG => CK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK2_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D0_CK2(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D0_CK2
!#endif
!        use pm_kind, only: CKG => CK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK1_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D0_CK1(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D0_CK1
!#endif
!        use pm_kind, only: CKG => CK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D0_RK5(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D0_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D0_RK4(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D0_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D0_RK3(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D0_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D0_RK2(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D0_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D0_RK1(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D0_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)                    :: sep
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
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
!    module subroutine setSplitFixCusComCusIns_D0_D0_SK5(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D0_D0_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)                    :: array
!        character(*,SKG)            , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    module subroutine setSplitFixCusComCusIns_D0_D0_SK4(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D0_D0_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)                    :: array
!        character(*,SKG)            , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    module subroutine setSplitFixCusComCusIns_D0_D0_SK3(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D0_D0_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)                    :: array
!        character(*,SKG)            , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    module subroutine setSplitFixCusComCusIns_D0_D0_SK2(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D0_D0_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)                    :: array
!        character(*,SKG)            , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    module subroutine setSplitFixCusComCusIns_D0_D0_SK1(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D0_D0_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)                    :: array
!        character(*,SKG)            , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if SK5_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D0_SK5(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D0_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D0_SK4(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D0_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D0_SK3(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D0_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D0_SK2(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D0_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D0_SK1(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D0_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D0_IK5(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D0_IK5
!#endif
!        use pm_kind, only: IKG => IK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D0_IK4(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D0_IK4
!#endif
!        use pm_kind, only: IKG => IK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D0_IK3(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D0_IK3
!#endif
!        use pm_kind, only: IKG => IK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D0_IK2(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D0_IK2
!#endif
!        use pm_kind, only: IKG => IK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D0_IK1(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D0_IK1
!#endif
!        use pm_kind, only: IKG => IK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if LK5_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D0_LK5(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D0_LK5
!#endif
!        use pm_kind, only: LKG => LK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK4_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D0_LK4(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D0_LK4
!#endif
!        use pm_kind, only: LKG => LK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK3_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D0_LK3(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D0_LK3
!#endif
!        use pm_kind, only: LKG => LK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK2_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D0_LK2(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D0_LK2
!#endif
!        use pm_kind, only: LKG => LK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK1_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D0_LK1(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D0_LK1
!#endif
!        use pm_kind, only: LKG => LK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if CK5_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D0_CK5(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D0_CK5
!#endif
!        use pm_kind, only: CKG => CK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK4_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D0_CK4(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D0_CK4
!#endif
!        use pm_kind, only: CKG => CK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK3_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D0_CK3(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D0_CK3
!#endif
!        use pm_kind, only: CKG => CK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK2_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D0_CK2(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D0_CK2
!#endif
!        use pm_kind, only: CKG => CK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK1_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D0_CK1(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D0_CK1
!#endif
!        use pm_kind, only: CKG => CK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D0_RK5(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D0_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D0_RK4(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D0_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D0_RK3(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D0_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D0_RK2(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D0_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D0_RK1(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D0_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)                    :: sep
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
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
!    PURE module subroutine setSplitFixDefComDefIns_D1_D1_SK5(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D1_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D1_SK4(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D1_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D1_SK3(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D1_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D1_SK2(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D1_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D1_SK1(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D1_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D1_IK5(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D1_IK5
!#endif
!        use pm_kind, only: IKG => IK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D1_IK4(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D1_IK4
!#endif
!        use pm_kind, only: IKG => IK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D1_IK3(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D1_IK3
!#endif
!        use pm_kind, only: IKG => IK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D1_IK2(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D1_IK2
!#endif
!        use pm_kind, only: IKG => IK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D1_IK1(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D1_IK1
!#endif
!        use pm_kind, only: IKG => IK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if LK5_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D1_LK5(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D1_LK5
!#endif
!        use pm_kind, only: LKG => LK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK4_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D1_LK4(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D1_LK4
!#endif
!        use pm_kind, only: LKG => LK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK3_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D1_LK3(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D1_LK3
!#endif
!        use pm_kind, only: LKG => LK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK2_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D1_LK2(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D1_LK2
!#endif
!        use pm_kind, only: LKG => LK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK1_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D1_LK1(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D1_LK1
!#endif
!        use pm_kind, only: LKG => LK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if CK5_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D1_CK5(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D1_CK5
!#endif
!        use pm_kind, only: CKG => CK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK4_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D1_CK4(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D1_CK4
!#endif
!        use pm_kind, only: CKG => CK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK3_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D1_CK3(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D1_CK3
!#endif
!        use pm_kind, only: CKG => CK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK2_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D1_CK2(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D1_CK2
!#endif
!        use pm_kind, only: CKG => CK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK1_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D1_CK1(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D1_CK1
!#endif
!        use pm_kind, only: CKG => CK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D1_RK5(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D1_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D1_RK4(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D1_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D1_RK3(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D1_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D1_RK2(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D1_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE module subroutine setSplitFixDefComDefIns_D1_D1_RK1(field, nsplit, array, sep, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComDefIns_D1_D1_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
!        logical(LK)                 , intent(in)    , optional      :: keep
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
!    module subroutine setSplitFixCusComDefIns_D1_D1_SK5(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D1_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D1_SK4(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D1_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D1_SK3(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D1_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D1_SK2(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D1_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D1_SK1(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D1_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D1_IK5(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D1_IK5
!#endif
!        use pm_kind, only: IKG => IK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D1_IK4(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D1_IK4
!#endif
!        use pm_kind, only: IKG => IK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D1_IK3(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D1_IK3
!#endif
!        use pm_kind, only: IKG => IK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D1_IK2(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D1_IK2
!#endif
!        use pm_kind, only: IKG => IK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D1_IK1(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D1_IK1
!#endif
!        use pm_kind, only: IKG => IK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if LK5_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D1_LK5(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D1_LK5
!#endif
!        use pm_kind, only: LKG => LK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK4_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D1_LK4(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D1_LK4
!#endif
!        use pm_kind, only: LKG => LK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK3_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D1_LK3(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D1_LK3
!#endif
!        use pm_kind, only: LKG => LK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK2_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D1_LK2(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D1_LK2
!#endif
!        use pm_kind, only: LKG => LK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK1_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D1_LK1(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D1_LK1
!#endif
!        use pm_kind, only: LKG => LK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if CK5_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D1_CK5(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D1_CK5
!#endif
!        use pm_kind, only: CKG => CK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK4_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D1_CK4(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D1_CK4
!#endif
!        use pm_kind, only: CKG => CK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK3_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D1_CK3(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D1_CK3
!#endif
!        use pm_kind, only: CKG => CK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK2_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D1_CK2(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D1_CK2
!#endif
!        use pm_kind, only: CKG => CK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK1_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D1_CK1(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D1_CK1
!#endif
!        use pm_kind, only: CKG => CK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D1_RK5(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D1_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D1_RK4(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D1_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D1_RK3(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D1_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D1_RK2(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D1_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    module subroutine setSplitFixCusComDefIns_D1_D1_RK1(field, nsplit, array, sep, iseq, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComDefIns_D1_D1_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        logical(LK)                 , intent(in)    , optional      :: keep
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
!    PURE module subroutine setSplitFixDefComCusIns_D1_D1_SK5(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D1_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D1_SK4(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D1_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D1_SK3(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D1_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D1_SK2(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D1_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D1_SK1(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D1_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D1_IK5(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D1_IK5
!#endif
!        use pm_kind, only: IKG => IK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D1_IK4(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D1_IK4
!#endif
!        use pm_kind, only: IKG => IK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D1_IK3(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D1_IK3
!#endif
!        use pm_kind, only: IKG => IK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D1_IK2(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D1_IK2
!#endif
!        use pm_kind, only: IKG => IK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D1_IK1(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D1_IK1
!#endif
!        use pm_kind, only: IKG => IK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if LK5_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D1_LK5(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D1_LK5
!#endif
!        use pm_kind, only: LKG => LK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK4_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D1_LK4(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D1_LK4
!#endif
!        use pm_kind, only: LKG => LK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK3_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D1_LK3(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D1_LK3
!#endif
!        use pm_kind, only: LKG => LK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK2_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D1_LK2(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D1_LK2
!#endif
!        use pm_kind, only: LKG => LK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK1_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D1_LK1(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D1_LK1
!#endif
!        use pm_kind, only: LKG => LK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if CK5_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D1_CK5(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D1_CK5
!#endif
!        use pm_kind, only: CKG => CK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK4_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D1_CK4(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D1_CK4
!#endif
!        use pm_kind, only: CKG => CK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK3_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D1_CK3(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D1_CK3
!#endif
!        use pm_kind, only: CKG => CK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK2_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D1_CK2(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D1_CK2
!#endif
!        use pm_kind, only: CKG => CK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK1_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D1_CK1(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D1_CK1
!#endif
!        use pm_kind, only: CKG => CK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D1_RK5(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D1_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D1_RK4(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D1_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D1_RK3(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D1_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D1_RK2(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D1_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE module subroutine setSplitFixDefComCusIns_D1_D1_RK1(field, nsplit, array, sep, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixDefComCusIns_D1_D1_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
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
!    module subroutine setSplitFixCusComCusIns_D1_D1_SK5(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D1_SK5
!#endif
!        use pm_kind, only: SKG => SK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK4_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D1_SK4(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D1_SK4
!#endif
!        use pm_kind, only: SKG => SK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK3_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D1_SK3(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D1_SK3
!#endif
!        use pm_kind, only: SKG => SK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK2_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D1_SK2(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D1_SK2
!#endif
!        use pm_kind, only: SKG => SK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if SK1_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D1_SK1(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D1_SK1
!#endif
!        use pm_kind, only: SKG => SK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
!        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK5_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D1_IK5(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D1_IK5
!#endif
!        use pm_kind, only: IKG => IK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK4_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D1_IK4(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D1_IK4
!#endif
!        use pm_kind, only: IKG => IK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK3_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D1_IK3(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D1_IK3
!#endif
!        use pm_kind, only: IKG => IK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK2_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D1_IK2(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D1_IK2
!#endif
!        use pm_kind, only: IKG => IK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if IK1_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D1_IK1(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D1_IK1
!#endif
!        use pm_kind, only: IKG => IK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        integer(IKG)                , intent(in)    , contiguous    :: array(:)
!        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if LK5_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D1_LK5(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D1_LK5
!#endif
!        use pm_kind, only: LKG => LK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK4_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D1_LK4(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D1_LK4
!#endif
!        use pm_kind, only: LKG => LK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK3_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D1_LK3(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D1_LK3
!#endif
!        use pm_kind, only: LKG => LK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK2_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D1_LK2(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D1_LK2
!#endif
!        use pm_kind, only: LKG => LK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if LK1_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D1_LK1(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D1_LK1
!#endif
!        use pm_kind, only: LKG => LK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        logical(LKG)                , intent(in)    , contiguous    :: array(:)
!        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if CK5_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D1_CK5(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D1_CK5
!#endif
!        use pm_kind, only: CKG => CK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK4_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D1_CK4(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D1_CK4
!#endif
!        use pm_kind, only: CKG => CK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK3_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D1_CK3(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D1_CK3
!#endif
!        use pm_kind, only: CKG => CK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK2_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D1_CK2(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D1_CK2
!#endif
!        use pm_kind, only: CKG => CK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if CK1_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D1_CK1(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D1_CK1
!#endif
!        use pm_kind, only: CKG => CK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        complex(CKG)                , intent(in)    , contiguous    :: array(:)
!        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D1_RK5(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D1_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D1_RK4(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D1_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D1_RK3(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D1_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D1_RK2(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D1_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    module subroutine setSplitFixCusComCusIns_D1_D1_RK1(field, nsplit, array, sep, iseq, instance, sorted, unique, keep)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitFixCusComCusIns_D1_D1_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        integer(IK)                 , intent(out)   , contiguous    :: field(:)
!        integer(IK)                 , intent(out)                   :: nsplit
!        real(RKG)                   , intent(in)    , contiguous    :: array(:)
!        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
!        procedure(logical(LK))                                      :: iseq
!        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
!        logical(LK)                 , intent(in)    , optional      :: sorted
!        logical(LK)                 , intent(in)    , optional      :: unique
!        logical(LK)                 , intent(in)    , optional      :: keep
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface

    ! allocatable index

    interface setSplit

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D0_D0_SK5(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D0_D0_SK4(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D0_D0_SK3(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D0_D0_SK2(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D0_D0_SK1(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D0_SK5(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D0_SK4(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D0_SK3(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D0_SK2(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D0_SK1(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D0_IK5(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D0_IK4(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D0_IK3(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D0_IK2(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D0_IK1(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D0_LK5(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D0_LK4(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D0_LK3(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D0_LK2(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D0_LK1(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D0_CK5(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D0_CK4(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D0_CK3(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D0_CK2(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D0_CK1(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D0_RK5(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D0_RK4(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D0_RK3(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D0_RK2(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D0_RK1(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSplitIndCusComDefIns_D0_D0_SK5(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSplitIndCusComDefIns_D0_D0_SK4(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSplitIndCusComDefIns_D0_D0_SK3(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSplitIndCusComDefIns_D0_D0_SK2(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSplitIndCusComDefIns_D0_D0_SK1(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D0_SK5(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D0_SK4(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D0_SK3(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D0_SK2(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D0_SK1(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D0_IK5(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D0_IK4(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D0_IK3(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D0_IK2(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D0_IK1(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D0_LK5(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D0_LK4(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D0_LK3(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D0_LK2(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D0_LK1(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D0_CK5(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D0_CK4(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D0_CK3(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D0_CK2(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D0_CK1(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D0_RK5(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D0_RK4(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D0_RK3(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D0_RK2(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D0_RK1(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D0_D0_SK5(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D0_D0_SK4(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D0_D0_SK3(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D0_D0_SK2(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D0_D0_SK1(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D0_SK5(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D0_SK4(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D0_SK3(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D0_SK2(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D0_SK1(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D0_IK5(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D0_IK4(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D0_IK3(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D0_IK2(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D0_IK1(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D0_LK5(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D0_LK4(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D0_LK3(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D0_LK2(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D0_LK1(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D0_CK5(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D0_CK4(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D0_CK3(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D0_CK2(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D0_CK1(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D0_RK5(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D0_RK4(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D0_RK3(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D0_RK2(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D0_RK1(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSplitIndCusComCusIns_D0_D0_SK5(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSplitIndCusComCusIns_D0_D0_SK4(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSplitIndCusComCusIns_D0_D0_SK3(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSplitIndCusComCusIns_D0_D0_SK2(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSplitIndCusComCusIns_D0_D0_SK1(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D0_SK5(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D0_SK4(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D0_SK3(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D0_SK2(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D0_SK1(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D0_IK5(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D0_IK4(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D0_IK3(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D0_IK2(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D0_IK1(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D0_LK5(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D0_LK4(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D0_LK3(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D0_LK2(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D0_LK1(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D0_CK5(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D0_CK4(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D0_CK3(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D0_CK2(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D0_CK1(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D0_RK5(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D0_RK4(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D0_RK3(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D0_RK2(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D0_RK1(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D1_SK5(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D1_SK4(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D1_SK3(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D1_SK2(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D1_SK1(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D1_IK5(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D1_IK4(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D1_IK3(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D1_IK2(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D1_IK1(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D1_LK5(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D1_LK4(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D1_LK3(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D1_LK2(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D1_LK1(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D1_CK5(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D1_CK4(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D1_CK3(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D1_CK2(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D1_CK1(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D1_RK5(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D1_RK4(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D1_RK3(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D1_RK2(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setSplitIndDefComDefIns_D1_D1_RK1(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComDefIns_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D1_SK5(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D1_SK4(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D1_SK3(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D1_SK2(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D1_SK1(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D1_IK5(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D1_IK4(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D1_IK3(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D1_IK2(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D1_IK1(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D1_LK5(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D1_LK4(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D1_LK3(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D1_LK2(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D1_LK1(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D1_CK5(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D1_CK4(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D1_CK3(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D1_CK2(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D1_CK1(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D1_RK5(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D1_RK4(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D1_RK3(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D1_RK2(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setSplitIndCusComDefIns_D1_D1_RK1(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComDefIns_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D1_SK5(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D1_SK4(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D1_SK3(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D1_SK2(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D1_SK1(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D1_IK5(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D1_IK4(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D1_IK3(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D1_IK2(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D1_IK1(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D1_LK5(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D1_LK4(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D1_LK3(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D1_LK2(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D1_LK1(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D1_CK5(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D1_CK4(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D1_CK3(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D1_CK2(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D1_CK1(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D1_RK5(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D1_RK4(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D1_RK3(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D1_RK2(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setSplitIndDefComCusIns_D1_D1_RK1(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndDefComCusIns_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D1_SK5(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D1_SK4(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D1_SK3(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D1_SK2(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D1_SK1(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D1_IK5(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D1_IK4(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D1_IK3(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D1_IK2(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D1_IK1(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D1_LK5(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D1_LK4(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D1_LK3(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D1_LK2(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D1_LK1(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D1_CK5(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D1_CK4(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D1_CK3(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D1_CK2(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D1_CK1(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D1_RK5(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D1_RK4(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D1_RK3(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D1_RK2(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setSplitIndCusComCusIns_D1_D1_RK1(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitIndCusComCusIns_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)                 , intent(out)   , allocatable   :: field(:,:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! PDT container

#if PDT_ENABLED
    interface setSplit

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D0_D0_SK5(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D0_D0_SK4(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D0_D0_SK3(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D0_D0_SK2(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D0_D0_SK1(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D0_SK5(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D0_SK4(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D0_SK3(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D0_SK2(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D0_SK1(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D0_IK5(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D0_IK4(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D0_IK3(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D0_IK2(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D0_IK1(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D0_LK5(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D0_LK4(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D0_LK3(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D0_LK2(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D0_LK1(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D0_CK5(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D0_CK4(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D0_CK3(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D0_CK2(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D0_CK1(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D0_RK5(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D0_RK4(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D0_RK3(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D0_RK2(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D0_RK1(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSplitConCusComDefIns_D0_D0_SK5(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSplitConCusComDefIns_D0_D0_SK4(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSplitConCusComDefIns_D0_D0_SK3(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSplitConCusComDefIns_D0_D0_SK2(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSplitConCusComDefIns_D0_D0_SK1(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D0_SK5(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D0_SK4(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D0_SK3(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D0_SK2(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D0_SK1(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D0_IK5(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D0_IK4(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D0_IK3(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D0_IK2(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D0_IK1(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D0_LK5(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D0_LK4(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D0_LK3(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D0_LK2(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D0_LK1(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D0_CK5(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D0_CK4(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D0_CK3(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D0_CK2(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D0_CK1(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D0_RK5(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D0_RK4(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D0_RK3(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D0_RK2(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D0_RK1(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D0_D0_SK5(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D0_D0_SK4(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D0_D0_SK3(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D0_D0_SK2(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D0_D0_SK1(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D0_SK5(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D0_SK4(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D0_SK3(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D0_SK2(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D0_SK1(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D0_IK5(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D0_IK4(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D0_IK3(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D0_IK2(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D0_IK1(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D0_LK5(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D0_LK4(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D0_LK3(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D0_LK2(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D0_LK1(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D0_CK5(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D0_CK4(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D0_CK3(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D0_CK2(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D0_CK1(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D0_RK5(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D0_RK4(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D0_RK3(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D0_RK2(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D0_RK1(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSplitConCusComCusIns_D0_D0_SK5(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSplitConCusComCusIns_D0_D0_SK4(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSplitConCusComCusIns_D0_D0_SK3(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSplitConCusComCusIns_D0_D0_SK2(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSplitConCusComCusIns_D0_D0_SK1(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D0_SK5(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D0_SK4(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D0_SK3(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D0_SK2(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D0_SK1(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D0_IK5(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D0_IK4(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D0_IK3(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D0_IK2(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D0_IK1(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D0_LK5(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D0_LK4(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D0_LK3(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D0_LK2(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D0_LK1(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D0_CK5(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D0_CK4(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D0_CK3(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D0_CK2(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D0_CK1(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D0_RK5(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D0_RK4(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D0_RK3(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D0_RK2(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D0_RK1(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D1_SK5(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D1_SK4(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D1_SK3(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D1_SK2(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D1_SK1(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D1_IK5(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D1_IK4(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D1_IK3(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D1_IK2(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D1_IK1(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D1_LK5(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D1_LK4(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D1_LK3(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D1_LK2(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D1_LK1(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D1_CK5(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D1_CK4(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D1_CK3(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D1_CK2(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D1_CK1(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D1_RK5(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D1_RK4(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D1_RK3(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D1_RK2(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setSplitConDefComDefIns_D1_D1_RK1(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComDefIns_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D1_SK5(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D1_SK4(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D1_SK3(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D1_SK2(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D1_SK1(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D1_IK5(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D1_IK4(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D1_IK3(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D1_IK2(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D1_IK1(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D1_LK5(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D1_LK4(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D1_LK3(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D1_LK2(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D1_LK1(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D1_CK5(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D1_CK4(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D1_CK3(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D1_CK2(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D1_CK1(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D1_RK5(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D1_RK4(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D1_RK3(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D1_RK2(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setSplitConCusComDefIns_D1_D1_RK1(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComDefIns_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D1_SK5(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D1_SK4(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D1_SK3(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D1_SK2(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D1_SK1(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D1_IK5(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D1_IK4(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D1_IK3(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D1_IK2(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D1_IK1(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D1_LK5(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D1_LK4(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D1_LK3(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D1_LK2(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D1_LK1(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D1_CK5(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D1_CK4(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D1_CK3(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D1_CK2(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D1_CK1(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D1_RK5(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D1_RK4(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D1_RK3(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D1_RK2(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setSplitConDefComCusIns_D1_D1_RK1(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConDefComCusIns_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D1_SK5(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D1_SK4(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D1_SK3(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D1_SK2(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D1_SK1(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: cvs_pdt
        type(cvs_pdt(SKG))          , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D1_IK5(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D1_IK4(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D1_IK3(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D1_IK2(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D1_IK1(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        use pm_container, only: cvi_pdt
        type(cvi_pdt(IKG))          , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D1_LK5(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D1_LK4(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D1_LK3(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D1_LK2(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D1_LK1(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        use pm_container, only: cvl_pdt
        type(cvl_pdt(LKG))          , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D1_CK5(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D1_CK4(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D1_CK3(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D1_CK2(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D1_CK1(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        use pm_container, only: cvc_pdt
        type(cvc_pdt(CKG))          , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D1_RK5(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D1_RK4(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D1_RK3(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D1_RK2(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setSplitConCusComCusIns_D1_D1_RK1(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitConCusComCusIns_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        use pm_container, only: cvr_pdt
        type(cvr_pdt(RKG))          , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface
#endif
!PDT_ENABLED

    ! Box container

    interface setSplit

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setSplitBoxDefComDefIns_D0_D0_SK(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxDefComDefIns_D0_D0_SK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setSplitBoxDefComDefIns_D1_D0_SK(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxDefComDefIns_D1_D0_SK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: cvs_type
        type(cvs_type)              , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setSplitBoxDefComDefIns_D1_D0_IK(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxDefComDefIns_D1_D0_IK
#endif
        use pm_kind, only: IKG => IK
        use pm_container, only: cvi_type
        type(cvi_type)              , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setSplitBoxDefComDefIns_D1_D0_LK(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxDefComDefIns_D1_D0_LK
#endif
        use pm_kind, only: LKG => LK
        use pm_container, only: cvl_type
        type(cvl_type)              , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setSplitBoxDefComDefIns_D1_D0_CK(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxDefComDefIns_D1_D0_CK
#endif
        use pm_kind, only: CKG => CK
        use pm_container, only: cvc_type
        type(cvc_type)              , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setSplitBoxDefComDefIns_D1_D0_RK(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxDefComDefIns_D1_D0_RK
#endif
        use pm_kind, only: RKG => RK
        use pm_container, only: cvr_type
        type(cvr_type)              , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSplitBoxCusComDefIns_D0_D0_SK(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxCusComDefIns_D0_D0_SK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSplitBoxCusComDefIns_D1_D0_SK(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxCusComDefIns_D1_D0_SK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: cvs_type
        type(cvs_type)              , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSplitBoxCusComDefIns_D1_D0_IK(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxCusComDefIns_D1_D0_IK
#endif
        use pm_kind, only: IKG => IK
        use pm_container, only: cvi_type
        type(cvi_type)              , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSplitBoxCusComDefIns_D1_D0_LK(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxCusComDefIns_D1_D0_LK
#endif
        use pm_kind, only: LKG => LK
        use pm_container, only: cvl_type
        type(cvl_type)              , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSplitBoxCusComDefIns_D1_D0_CK(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxCusComDefIns_D1_D0_CK
#endif
        use pm_kind, only: CKG => CK
        use pm_container, only: cvc_type
        type(cvc_type)              , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSplitBoxCusComDefIns_D1_D0_RK(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxCusComDefIns_D1_D0_RK
#endif
        use pm_kind, only: RKG => RK
        use pm_container, only: cvr_type
        type(cvr_type)              , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setSplitBoxDefComCusIns_D0_D0_SK(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxDefComCusIns_D0_D0_SK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setSplitBoxDefComCusIns_D1_D0_SK(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxDefComCusIns_D1_D0_SK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: cvs_type
        type(cvs_type)              , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setSplitBoxDefComCusIns_D1_D0_IK(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxDefComCusIns_D1_D0_IK
#endif
        use pm_kind, only: IKG => IK
        use pm_container, only: cvi_type
        type(cvi_type)              , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setSplitBoxDefComCusIns_D1_D0_LK(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxDefComCusIns_D1_D0_LK
#endif
        use pm_kind, only: LKG => LK
        use pm_container, only: cvl_type
        type(cvl_type)              , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setSplitBoxDefComCusIns_D1_D0_CK(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxDefComCusIns_D1_D0_CK
#endif
        use pm_kind, only: CKG => CK
        use pm_container, only: cvc_type
        type(cvc_type)              , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setSplitBoxDefComCusIns_D1_D0_RK(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxDefComCusIns_D1_D0_RK
#endif
        use pm_kind, only: RKG => RK
        use pm_container, only: cvr_type
        type(cvr_type)              , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSplitBoxCusComCusIns_D0_D0_SK(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxCusComCusIns_D0_D0_SK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)                    :: array
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSplitBoxCusComCusIns_D1_D0_SK(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxCusComCusIns_D1_D0_SK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: cvs_type
        type(cvs_type)              , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSplitBoxCusComCusIns_D1_D0_IK(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxCusComCusIns_D1_D0_IK
#endif
        use pm_kind, only: IKG => IK
        use pm_container, only: cvi_type
        type(cvi_type)              , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSplitBoxCusComCusIns_D1_D0_LK(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxCusComCusIns_D1_D0_LK
#endif
        use pm_kind, only: LKG => LK
        use pm_container, only: cvl_type
        type(cvl_type)              , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSplitBoxCusComCusIns_D1_D0_CK(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxCusComCusIns_D1_D0_CK
#endif
        use pm_kind, only: CKG => CK
        use pm_container, only: cvc_type
        type(cvc_type)              , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSplitBoxCusComCusIns_D1_D0_RK(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxCusComCusIns_D1_D0_RK
#endif
        use pm_kind, only: RKG => RK
        use pm_container, only: cvr_type
        type(cvr_type)              , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)                    :: sep
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setSplitBoxDefComDefIns_D1_D1_SK(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxDefComDefIns_D1_D1_SK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: cvs_type
        type(cvs_type)              , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setSplitBoxDefComDefIns_D1_D1_IK(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxDefComDefIns_D1_D1_IK
#endif
        use pm_kind, only: IKG => IK
        use pm_container, only: cvi_type
        type(cvi_type)              , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setSplitBoxDefComDefIns_D1_D1_LK(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxDefComDefIns_D1_D1_LK
#endif
        use pm_kind, only: LKG => LK
        use pm_container, only: cvl_type
        type(cvl_type)              , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setSplitBoxDefComDefIns_D1_D1_CK(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxDefComDefIns_D1_D1_CK
#endif
        use pm_kind, only: CKG => CK
        use pm_container, only: cvc_type
        type(cvc_type)              , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setSplitBoxDefComDefIns_D1_D1_RK(field, array, sep, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxDefComDefIns_D1_D1_RK
#endif
        use pm_kind, only: RKG => RK
        use pm_container, only: cvr_type
        type(cvr_type)              , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSplitBoxCusComDefIns_D1_D1_SK(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxCusComDefIns_D1_D1_SK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: cvs_type
        type(cvs_type)              , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSplitBoxCusComDefIns_D1_D1_IK(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxCusComDefIns_D1_D1_IK
#endif
        use pm_kind, only: IKG => IK
        use pm_container, only: cvi_type
        type(cvi_type)              , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSplitBoxCusComDefIns_D1_D1_LK(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxCusComDefIns_D1_D1_LK
#endif
        use pm_kind, only: LKG => LK
        use pm_container, only: cvl_type
        type(cvl_type)              , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSplitBoxCusComDefIns_D1_D1_CK(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxCusComDefIns_D1_D1_CK
#endif
        use pm_kind, only: CKG => CK
        use pm_container, only: cvc_type
        type(cvc_type)              , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSplitBoxCusComDefIns_D1_D1_RK(field, array, sep, iseq, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxCusComDefIns_D1_D1_RK
#endif
        use pm_kind, only: RKG => RK
        use pm_container, only: cvr_type
        type(cvr_type)              , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setSplitBoxDefComCusIns_D1_D1_SK(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxDefComCusIns_D1_D1_SK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: cvs_type
        type(cvs_type)              , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setSplitBoxDefComCusIns_D1_D1_IK(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxDefComCusIns_D1_D1_IK
#endif
        use pm_kind, only: IKG => IK
        use pm_container, only: cvi_type
        type(cvi_type)              , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setSplitBoxDefComCusIns_D1_D1_LK(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxDefComCusIns_D1_D1_LK
#endif
        use pm_kind, only: LKG => LK
        use pm_container, only: cvl_type
        type(cvl_type)              , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setSplitBoxDefComCusIns_D1_D1_CK(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxDefComCusIns_D1_D1_CK
#endif
        use pm_kind, only: CKG => CK
        use pm_container, only: cvc_type
        type(cvc_type)              , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setSplitBoxDefComCusIns_D1_D1_RK(field, array, sep, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxDefComCusIns_D1_D1_RK
#endif
        use pm_kind, only: RKG => RK
        use pm_container, only: cvr_type
        type(cvr_type)              , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSplitBoxCusComCusIns_D1_D1_SK(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxCusComCusIns_D1_D1_SK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: cvs_type
        type(cvs_type)              , intent(out)   , allocatable   :: field(:)
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        character(*,SKG)            , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSplitBoxCusComCusIns_D1_D1_IK(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxCusComCusIns_D1_D1_IK
#endif
        use pm_kind, only: IKG => IK
        use pm_container, only: cvi_type
        type(cvi_type)              , intent(out)   , allocatable   :: field(:)
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSplitBoxCusComCusIns_D1_D1_LK(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxCusComCusIns_D1_D1_LK
#endif
        use pm_kind, only: LKG => LK
        use pm_container, only: cvl_type
        type(cvl_type)              , intent(out)   , allocatable   :: field(:)
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        logical(LKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSplitBoxCusComCusIns_D1_D1_CK(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxCusComCusIns_D1_D1_CK
#endif
        use pm_kind, only: CKG => CK
        use pm_container, only: cvc_type
        type(cvc_type)              , intent(out)   , allocatable   :: field(:)
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        complex(CKG)                , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSplitBoxCusComCusIns_D1_D1_RK(field, array, sep, iseq, instance, sorted, unique, keep)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSplitBoxCusComCusIns_D1_D1_RK
#endif
        use pm_kind, only: RKG => RK
        use pm_container, only: cvr_type
        type(cvr_type)              , intent(out)   , allocatable   :: field(:)
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        real(RKG)                   , intent(in)    , contiguous    :: sep(:)
        procedure(logical(LK))                                      :: iseq
        integer(IK)                 , intent(in)    , contiguous    :: instance(:)
        logical(LK)                 , intent(in)    , optional      :: sorted
        logical(LK)                 , intent(in)    , optional      :: unique
        logical(LK)                 , intent(in)    , optional      :: keep
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arraySplit ! LCOV_EXCL_LINE