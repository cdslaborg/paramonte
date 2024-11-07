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
!>  This module contains a collection of example functions
!>  for testing or examining the root-finding routines of the ParaMonte library.
!>
!>  \see
!>  [pm_mathRoot](@ref pm_mathRoot)<br>
!>  [pm_mathRoot](@ref pm_mathRoot)<br>
!>  [pm_mathRoot](@ref pm_mathRoot)<br>
!>
!>  \todo
!>  \pvhigh
!>  Unfortunately, gfortran 12 and older versions do not properly support the parameterized derived types (PDTs).<br>
!>  As such, the example generic-real-kind PDT types could not be used here.<br>
!>  This creates significant complexities when using the examples of these modules,<br>
!>  because all `real` kinds in this module are set to the highest precision available.<br>
!>  The onus is then on the end user to write wrappers that convert the relevant components and function-returns to the desired `real` kinds.<br>
!>  In addition to being ugly, error-prone and verbose, this usage of the highest-precision `real` kind is also highly inefficient computationally.<br>
!>  Fortunately, once PDTs are supported in gfortran, the conversion of the example types of this module to PDTs is straightforward and non-breaking.<br>
!>  The migration to PDTs must be done as soon as gfortran supports for PDTs is complete.<br>
!>  Note that other compilers have full support of PDTs.<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_mathRootTest

    use pm_val2str, only: getStr
    use pm_option, only: getOption
    use pm_kind, only: SK, IK, LK, RKH
    use pm_str, only: getTTZ => getTrimmedTZ

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathRootTest"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the base type `func_type` standing `abstract` function type to generate a variety of function test functions.<br>
    !>
    !>  \details
    !>  The abstract type minimally contains the function limits and deferred type-bound procedure `get()`.<br>
    !>  It is primarily meant to be used internally within the ParaMonte library for testing and illustration purposes.<br>
    !>  In particular, the implementations may not be ideal for benchmarks and performance tests of the library integrators.<br>
    !>  Nevertheless, all relevant types and function test objects are `public` in this repository and available to the end user.<br>
    !>
    !>  \interface{func_type}
    !>  \code{.F90}
    !>
    !>      use pm_mathRootTest, only: func_type
    !>
    !>      type, extends(func_type) :: Func_type
    !>          ...
    !>      end type
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [getQuadErr](@ref pm_quadPack::getQuadErr)<br>
    !>
    !>  \final{func_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, abstract :: func_type
        real(RKH)                       :: lb       !<  \public The scalar of type `real` of the highest kind supported by the processor \RKH, containing the lower limit of function.
        real(RKH)                       :: ub       !<  \public The scalar of type `real` of the highest kind supported by the processor \RKH, containing the upper limit of function.
        real(RKH)       , allocatable   :: root(:)  !<  \public The scalar of type `real` of the highest kind supported by the processor \RKH, containing the true roots of function.
        character(:, SK), allocatable   :: desc     !<  \public The scalar `allocatable` character of default kind \SK containing a description of the function and function limits and difficulties.
    contains
        procedure(get_proc) , deferred  :: get
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the abstract interface of the `get()` type-bound procedure of [func_type](@ref pm_mathRootTest::func_type)
    !>  class whose arguments of type `real` are of the highest precision kind \RKH, made available by the processor.
    !>
    !>  \see
    !>  [func_type](@ref pm_mathRootTest::func_type)<br>
    !>
    !>  \final{get_proc}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    abstract interface
        PURE function get_proc(self, x) result(func)
            use pm_kind, only: RKG => RKH
            import :: func_type
            class(func_type)    , intent(in)    :: self
            real(RKG)           , intent(in)    :: x
            real(RKG)                           :: func
        end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for generating test function objects of the algebraic form as described below.
    !>
    !>  \details
    !>  The full function is defined as,<br>
    !>  \f{equation}{
    !>
    !>      f(x) = x (x^2 - 1)(x^2 - 4) ~,~ x \in (\ms{lb}, \ms{ub})
    !>
    !>  \f}
    !>  where the function bounds could be infinities.<br>
    !>  The real roots of the function are `-2, -1, 0, 1, 2`.<br>
    !>
    !>  \param[in]  lb      :   The input scalar of type `real` of kind \RKH, containing the lower limit of the root search bracket.<br>
    !>                          (**optional**, default = `-3.`)<br>
    !>  \param[in]  ub      :   The input scalar of the same type and kind as `lb`, containing the upper limit of the root search bracket.<br>
    !>                          (**optional**, default = `+3.`)<br>
    !>
    !>  \interface{func1_type}
    !>  \code{.F90}
    !>
    !>      use pm_mathRootTest, only: func1_type
    !>      type(func1_type) :: Func
    !>
    !>      Func = func1_type(lb = lb, ub = ub)
    !>      print *, "description: ", Func%desc
    !>      print *, "lower search limit: ", Func%lb
    !>      print *, "upper search limit: ", Func%ub
    !>      print *, "Example function value: ", Func%get(x)
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [func_type](@ref pm_mathRootTest::func_type)<br>
    !>
    !>  \test
    !>  [test_pm_quadPack](@ref test_pm_quadPack)
    !>
    !>  \final{func1_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    type, extends(func_type) :: func1_type
    contains
        procedure               :: get => getFunc1
    end type

    !>  \cond excluded
    interface func1_type
    PURE module function func1_typer(lb, ub) result(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: func1_typer
#endif
        use pm_kind, only: RKG => RKH
        real(RKG), intent(in), optional :: lb, ub
        type(func1_type)                :: self
    end function
    end interface
    !>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
    PURE module function getFunc1(self, x) result(func)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFunc1
#endif
        use pm_kind, only: RKG => RKH
        class(func1_type)   , intent(in)    :: self
        real(RKG)           , intent(in)    :: x
        real(RKG)                           :: func
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_mathRootTest