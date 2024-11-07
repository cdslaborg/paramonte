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
!>  This module contains a collection of interesting or challenging integrands
!>  for testing or examining the integration routines of the ParaMonte library.<br>
!>
!>  \details
!>  The routines to be tested include but are not limited to those of<br>
!>  <ul>
!>      <li>    [pm_quadPack](@ref pm_quadPack)<br>
!>      <li>    [pm_quadRomb](@ref pm_quadRomb)<br>
!>  </ul>
!>
!>  All test integrands are wrapped in a derived type of the base
!>  abstract class [integrand_type](@ref pm_quadTest::integrand_type).<br>
!>
!>  \see
!>  [pm_quadPack](@ref pm_quadPack)<br>
!>  [pm_quadRomb](@ref pm_quadRomb)<br>
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

module pm_quadTest

    use pm_kind, only: SK, IK, LK, RKH
    use pm_str, only: getTTZ => getTrimmedTZ
    use pm_quadPack, only: wcauchy_type
    use pm_option, only: getOption
    use pm_val2str, only: getStr

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_quadTest"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the base type `integrand_type` standing `abstract` integrand type to generate a variety of integration test functions.<br>
    !>
    !>  \details
    !>  The abstract type minimally contains:<br>
    !>  <ol>
    !>      <li>    The integration limits and deferred type-bound procedure `get()`.<br>
    !>      <li>    The locations of points of difficulties within the domain of integration
    !>              (other than those that appear in function weights such as Cauchy types of singularities).<br>
    !>      <li>    Information about the Cauchy type singularity and Cauchy weight of the function.<br>
    !>              The Cauchy weight has the form,<br>
    !>              \f{equation}{
    !>                  w(x) = \frac{1}{x - c}
    !>              \f}
    !>  </ol>
    !>  It is primarily meant to be used internally within the ParaMonte library for testing and illustration purposes.<br>
    !>  In particular, the implementations may not be ideal for benchmarks and performance tests of the library integrators.<br>
    !>  Nevertheless, all relevant types and integration test objects are `public` in this repository and available to the end user.<br>
    !>
    !>  \interface{integrand_type}
    !>  \code{.F90}
    !>
    !>      use pm_quadTest, only: integrand_type
    !>
    !>      type, extends(integrand_type) :: myint_type
    !>          ...
    !>      end type
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [getQuadErr](@ref pm_quadPack::getQuadErr)<br>
    !>
    !>  \test
    !>  [test_pm_quadPack](@ref test_pm_quadPack)
    !>
    !>  \final{integrand_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, abstract :: integrand_type
        real(RKH)                           :: lb           !<  \public The scalar of type `real` of the highest kind supported by the library \RKH, containing the lower limit of integration.
        real(RKH)                           :: ub           !<  \public The scalar of type `real` of the highest kind supported by the library \RKH, containing the upper limit of integration.
        real(RKH)                           :: integral     !<  \public The scalar of type `real` of the highest kind supported by the library \RKH, containing the true result of integration.
        real(RKH)           , allocatable   :: break(:)     !<  \public The scalar of type `real` of the highest kind supported by the library \RKH, containing the points of difficulties of integration.
        type(wcauchy_type)  , allocatable   :: wcauchy      !<  \public The scalar of type [wcauchy_type](@ref pm_quadPack::wcauchy_type), containing the Cauchy singularity of the integrand.
        character(:, SK)    , allocatable   :: desc         !<  \public The scalar `allocatable` character of default kind \SK containing a description of the integrand and integration limits and difficulties.
    contains
        procedure(get_proc) , deferred      :: get          !<  \public The function member returning the value of the <b>un</b>weighted integrand (whether Cauchy/sin/cos/algebraically types of weights) at a specified input point `x`.
    end type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the abstract interface of the `get()` type-bound procedure of [integrand_type](@ref pm_quadTest::integrand_type)
    !>  class whose arguments of type `real` are of the highest precision kind \RKH, made available by the processor.
    !>
    !>  \param[in]  x           :   The input scalar `real` of kind \RKH, containing the point at which the integrand must be computed.
    !   \param[in]  weighted    :   The input scalar `logical` of default kind \LK.
    !                               <ol>
    !                                   <li>    If it is `.true.`, the function value will be computed **including** its (Cauchy, sin, cos, algebraic, or other type of) weight.<br>
    !                                           Use this if the goal is to test the performance of the algorithms in handling points of difficulties without explicitly specifying them for the integrators.<br>
    !                                   <li>    If it is `.false.`, the function value will be computed **excluding** its (Cauchy, sin, cos, algebraic, or other type of) weight.<br>
    !                                           This is typically the value that should be passed to the integrators of [pm_quadPack](@ref pm_quadPack) module.<br>
    !                               </ol>
    !>
    !>  \see
    !>  [integrand_type](@ref pm_quadTest::integrand_type)<br>
    !>
    !>  \final{get_proc}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    abstract interface
        function get_proc(self, x) result(func)
            use pm_kind, only: RKG => RKH
            import :: integrand_type
            class(integrand_type)   , intent(in)    :: self
            real(RKG)               , intent(in)    :: x
            real(RKG)                               :: func
        end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Run the adaptive global quadrature methods for the specified input integrand object.
    !>
    !>  \details
    !>  This procedure is created solely for the purpose of facilitating the display of
    !>  the results of the quadrature of example integrands of [pm_quadTest](@ref pm_quadTest).
    !>
    !>
    !>  \param[inout]   disp        :   The input/output object of type [display_type](@ref pm_io::display_type) containing
    !>                                  information about the display on which the integration procedure and results should be displayed.<br>
    !>  \param[in]      integrand   :   The input object of class [integrand_type](@ref pm_quadTest::integrand_type)
    !>                                  containing information about the example integrand that is to be integrated.<br>
    !>  \param[in]      abstol      :   The input positive scalar of the type `real` of kind corresponding to the highest precision available \RKH,
    !>                                  representing the absolute tolerance for the integration result.<br>
    !>                                  (**optional**. The default value is set by [isFailedQuad]@(ref pm_quadPack::isFailedQuad))<br>
    !>  \param[in]      reltol      :   The input positive scalar of the type `real` of kind corresponding to the highest precision available \RKH,
    !>                                  representing the relative tolerance for the integration result.<br>
    !>                                  (**optional**, The default value is set by [isFailedQuad]@(ref pm_quadPack::isFailedQuad))<br>
    !>  \param[in]      nintmax     :   The input positive scalar of the type `integer` of default kind \IK,
    !>                                  representing the maximum number of adaptive interval formations allowed.<br>
    !>                                  (**optional**, default = `2000`)<br>
    !>
    !>  \interface{test_isFailedQuad}
    !>  \code{.F90}
    !>
    !>      use pm_quadTest, only: test_isFailedQuad
    !>
    !>      call test_isFailedQuad(disp, integrand, abstol = abstol, reltol = reltol, nintmax = nintmax)
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [integrand_type](@ref pm_quadTest::integrand_type)<br>
    !>
    !>  \test
    !>  [test_pm_quadPack](@ref test_pm_quadPack)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{11}
    !>  \desc
    !>  The \gfortran cannot handle **submodule** procedures with implicit procedures whose interfaces are **solely** declared in the parent module.<br>
    !>  The \gfortran cannot recognize the procedure arguments without duplicating the full interface in the submodule.<br>
    !>  For example, gfortran fails to compile the following submodule procedure interface,<br>
    !>  For example, gfortran returns the following error code,
    !>  \code{.sh}
    !>
    !>      75 |         call disp%show("integrand%desc")
    !>      |                          1
    !>      Error: ‘show’ at (1) should be a SUBROUTINE
    !>  \endcode
    !>  \remedy
    !>  Avoid this coding style until the bug is resolved.
    !>
    !>  \final{test_isFailedQuad}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface test_isFailedQuad
    module subroutine test_isFailedQuad_RKH(disp, integrand, abstol, reltol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: test_isFailedQuad_RKH
#endif
        use pm_kind, only: RKG => RKH
        use pm_io, only: display_type
        type(display_type)      , intent(inout)         :: disp
        class(integrand_type)   , intent(in)            :: integrand
        real(RKG)               , intent(in), optional  :: abstol, reltol
    end subroutine
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Run the adaptive global quadrature methods for the specified input integrand object.
    !>
    !>  \details
    !>  This procedure is created solely for the purpose of facilitating the display of
    !>  the results of the quadrature of example integrands of [pm_quadTest](@ref pm_quadTest).
    !>
    !>
    !>  \param[inout]   disp        :   The input/output object of type [display_type](@ref pm_io::display_type) containing
    !>                                  information about the display on which the integration procedure and results should be displayed.<br>
    !>  \param[in]      integrand   :   The input object of class [integrand_type](@ref pm_quadTest::integrand_type)
    !>                                  containing information about the example integrand that is to be integrated.<br>
    !>  \param[in]      atol        :   The input positive scalar of the type `real` of kind corresponding to the highest precision available \RKH,
    !>                                  representing the absolute tolerance for the integration result.<br>
    !>                                  (**optional**, default = `0._RKH`)<br>
    !>  \param[in]      rtol        :   The input positive scalar of the type `real` of kind corresponding to the highest precision available \RKH,
    !>                                  representing the relative tolerance for the integration result.<br>
    !>                                  (**optional**, default = `epsilon(0._RKH)**(2/3.)`)<br>
    !>  \param[in]      nintmax     :   The input positive scalar of the type `integer` of default kind \IK,
    !>                                  representing the maximum number of adaptive interval formations allowed.<br>
    !>                                  (**optional**, default = `2000`)<br>
    !>
    !>  \interface{test_getQuadErr}
    !>  \code{.F90}
    !>
    !>      use pm_quadTest, only: test_getQuadErr
    !>
    !>      call test_getQuadErr(disp, integrand, atol = atol, rtol = rtol, nintmax = nintmax)
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [integrand_type](@ref pm_quadTest::integrand_type)<br>
    !>
    !>  \test
    !>  [test_pm_quadPack](@ref test_pm_quadPack)
    !>
    !>  \final{test_getQuadErr}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface test_getQuadErr
    module subroutine test_getQuadErr_RKH(disp, integrand, atol, rtol, nintmax)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: test_getQuadErr_RKH
#endif
        use pm_kind, only: RKG => RKH
        use pm_io, only: display_type
        type(display_type)      , intent(inout)         :: disp
        class(integrand_type)   , intent(in)            :: integrand
        real(RKG)               , intent(in), optional  :: atol, rtol
        integer(IK)             , intent(in), optional  :: nintmax
    end subroutine
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for generating test integrand objects of the algebraic form as described below.
    !>
    !>  \details
    !>  The full integrand is defined as,<br>
    !>  \f{equation}{
    !>
    !>      f(x) = \frac{x^2}{(x^2 + 1)(x^2 + 4)} ~,~ x \in (\ms{lb}, \ms{ub})
    !>
    !>  \f}
    !>  where the integration bounds could be infinities.<br>
    !>  The indefinite integral of the integrand is,<br>
    !>  \f{equation}{
    !>
    !>      \int f(x) dx = \frac{2}{3} \tan^{-1}\bigg(\frac{x}{2}\bigg) - \frac{1}{3} \tan^{-1}(x) + \mathrm{constant} ~.
    !>
    !>  \f}
    !>
    !>
    !>  \param[in]  lb      :   The input scalar of type `real` of kind \RKH, containing the lower limit of integration.<br>
    !>                          (**optional**, default = [getInfNeg(real(0,kind(lb))](@ref pm_except::getInfNeg))<br>
    !>  \param[in]  ub      :   The input scalar of the same type and kind as `lb`, containing the upper limit of integration.<br>
    !>                          (**optional**, default = [getInfPos(real(0,kind(ub))](@ref pm_except::getInfPos))<br>
    !>
    !>  \interface{int1_type}
    !>  \code{.F90}
    !>
    !>      use pm_quadTest, only: int1_type
    !>      type(int1_type) :: integrand
    !>
    !>      integrand = int1_type(lb = lb, ub = ub)
    !>      print *, "description: ", integrand%desc
    !>      print *, "lower limit: ", integrand%lb
    !>      print *, "upper limit: ", integrand%ub
    !>      print *, "Example integrand value: ", integrand%get(x)
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [integrand_type](@ref pm_quadTest::integrand_type)<br>
    !>
    !>  \test
    !>  [test_pm_quadPack](@ref test_pm_quadPack)
    !>
    !>  \final{int1_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    type, extends(integrand_type)   :: int1_type
    contains
        procedure                   :: get => getInt1
    end type

    !>  \cond excluded
    interface int1_type
    module function int1_typer(lb, ub) result(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: int1_typer
#endif
        use pm_kind, only: RKG => RKH
        real(RKG), intent(in), optional :: lb, ub
        type(int1_type)                 :: self
    end function
    end interface
    !>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface
    module function getInt1(self, x) result(func)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInt1
#endif
        use pm_kind, only: RKG => RKH
        class(int1_type)    , intent(in)    :: self
        real(RKG)           , intent(in)    :: x
        real(RKG)                           :: func
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for generating test integrand objects of algebraic form as described below.
    !>
    !>  \details
    !>  The full integrand is defined over a finite interval as,<br>
    !>  \f{equation}{
    !>
    !>      f(x) = \frac{1}{\sqrt{a - bx}} ~,~ x \in (0, a / b) ~,~ a > 0 ~,~ b > 0 ~,
    !>
    !>  \f}
    !>  where the factors \f$a\f$ and \f$b\f$ are any finite positive real numbers.<br>
    !>  The integrand has a singularity at the upper bound of integration and has the precise value,<br>
    !>  \f{equation}{
    !>
    !>      \int_{\ms{lb} = 0}^{\ms{ub} = a / b} f(x) dx = -\frac{(2 \sqrt{a - bx})}{b} \bigg|_{lb}^{ub} ~.
    !>
    !>  \f}
    !>
    !>  \param[in]  a       :   The input positive-valued scalar of type `real` of kind \RKH, such that `a / b` represents the upper limit of integration.<br>
    !>                          (**optional**, default = `1.`)<br>
    !>  \param[in]  b       :   The input positive-valued scalar of type `real` of kind \RKH, such that `a / b` represents the upper limit of integration.<br>
    !>                          (**optional**, default = `1.`)<br>
    !>
    !>  \interface{int2_type}
    !>  \code{.F90}
    !>
    !>      use pm_quadTest, only: int2_type
    !>      type(intSinCos_type) :: integrand
    !>
    !>      integrand = int2_type(a = a, b = b)
    !>      print *, "description: ", integrand%desc
    !>      print *, "lower limit: ", integrand%lb
    !>      print *, "upper limit: ", integrand%ub
    !>      print *, "Example integrand value: ", integrand%get(x)
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [integrand_type](@ref pm_quadTest::integrand_type)<br>
    !>
    !>  \test
    !>  [test_pm_quadPack](@ref test_pm_quadPack)
    !>
    !>  \final{int2_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    type, extends(integrand_type)   :: int2_type
        real(RKH)                   :: a, b
    contains
        procedure                   :: get => getInt2
    end type

    !>  \cond excluded
    interface int2_type
    module function int2_typer(a, b) result(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: int2_typer
#endif
        use pm_kind, only: RKG => RKH
        real(RKG), intent(in), optional :: a, b
        type(int2_type)                 :: self
    end function
    end interface
    !>  \endcond excluded

    interface
    module function getInt2(self, x) result(func)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInt2
#endif
        use pm_kind, only: RKG => RKH
        class(int2_type)    , intent(in)    :: self
        real(RKG)           , intent(in)    :: x
        real(RKG)                           :: func
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for generating test integrand objects of algebraic form as described below.
    !>
    !>  \details
    !>  The full integrand is defined over a finite interval as,<br>
    !>  \f{equation}{
    !>
    !>      f(x) = \frac{\log(x)}{\sqrt{x}} ~,~ x \in (\ms{lb} = 0, \ms{ub}) ~.
    !>
    !>  \f}
    !>  The integrand has the precise value,<br>
    !>  \f{equation}{
    !>
    !>      \int_{\ms{lb} = 0}^{\ms{ub}} f(x) dx = 2 \sqrt{\ms{ub}} (\log(\ms{ub}) - 2) ~.
    !>
    !>  \f}
    !>
    !>  \param[in]  ub      :   The input positive scalar of type `real` of kind \RKH, containing the upper limit of integration.<br>
    !>                          (**optional**, default = `1.`)<br>
    !>
    !>  \interface{int3_type}
    !>  \code{.F90}
    !>
    !>      use pm_quadTest, only: int3_type
    !>      type(intSinCos_type) :: integrand
    !>
    !>      integrand = int3_type(ub = ub)
    !>      print *, "description: ", integrand%desc
    !>      print *, "lower limit: ", integrand%lb
    !>      print *, "upper limit: ", integrand%ub
    !>      print *, "Example integrand value: ", integrand%get(x)
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [integrand_type](@ref pm_quadTest::integrand_type)<br>
    !>
    !>  \test
    !>  [test_pm_quadPack](@ref test_pm_quadPack)
    !>
    !>  \final{int3_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    type, extends(integrand_type)   :: int3_type
    contains
        procedure                   :: get => getInt3
    end type

    !>  \cond excluded
    interface int3_type
    module function int3_typer(ub) result(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: int3_typer
#endif
        use pm_kind, only: RKG => RKH
        real(RKG), intent(in), optional :: ub
        type(int3_type)                 :: self
    end function
    end interface
    !>  \endcond excluded

    interface
    module function getInt3(self, x) result(func)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInt3
#endif
        use pm_kind, only: RKG => RKH
        class(int3_type)    , intent(in)    :: self
        real(RKG)           , intent(in)    :: x
        real(RKG)                           :: func
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for generating test integrand objects of the following algebraic form.
    !>
    !>  \details
    !>  The full integrand is defined as,
    !>  \f{equation}{
    !>
    !>      f(x) = \frac{\log(x)}{(\log(x)^2 + 1)^2} , x \in (0, 1)
    !>
    !>  \f}
    !>  with an integral of `-0.189275187882093321180367135892330338053417661540147291526012234`.<br>
    !>  This integrand is inspired by the examples of John Burkardt test suite for QAGWS routine of QuadPack.<br>
    !>
    !>  \interface{int4_type}
    !>  \code{.F90}
    !>
    !>      use pm_quadTest, only: int4_type
    !>      type(int4_type) :: integrand
    !>
    !>      integrand = int4_type()
    !>      print *, "description: ", integrand%desc
    !>      print *, "lower limit: ", integrand%lb
    !>      print *, "upper limit: ", integrand%ub
    !>      print *, "Example integrand value: ", integrand%get(x)
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [integrand_type](@ref pm_quadTest::integrand_type)<br>
    !>
    !>  \test
    !>  [test_pm_quadPack](@ref test_pm_quadPack)
    !>
    !>  \final{int4_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    type, extends(integrand_type)   :: int4_type
    contains
        procedure                   :: get => getInt4
    end type

    !>  \cond excluded
    interface int4_type
    module function int4_typer() result(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: int4_typer
#endif
        use pm_kind, only: RKG => RKH
        type(int4_type) :: self
    end function
    end interface
    !>  \endcond excluded

    interface
    module function getInt4(self, x) result(func)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInt4
#endif
        use pm_kind, only: RKG => RKH
        class(int4_type), intent(in)    :: self
        real(RKG)       , intent(in)    :: x
        real(RKG)                       :: func
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for generating test integrand objects of the following algebraic form.
    !>
    !>  \details
    !>  The full integrand is defined as,
    !>  \f{equation}{
    !>
    !>      f(x) = x^3 \log(\big|(x^2 - 1) (x^2 - 2.)\big|) ~,~ x \in (\ms{lb}, \ms{ub})
    !>
    !>  \f}
    !>  with four possible singularities depending on the choice of integration range: \f$ [ -\sqrt{2}, -1, 1, \sqrt{2} ] \f$<br>
    !>  The integral is of the form,<br>
    !>  \f{equation}{
    !>
    !>      \int_{\ms{lb}}^{\ms{lb}} f(x) dx = 0.25 \bigg[ x^4 \log( \big|(x^2 - 1) (x^2 - 2)\big| ) ) - 4\log(x^2 - 2) - \log(x^2 - 1) - 3x^2 - x^4 \bigg] \bigg|_{\ms{lb}}^{\ms{ub}}
    !>
    !>  \f}
    !>  This integrand is inspired by and extends the examples of John Burkardt test suite for QAGWS routine of QuadPack.<br>
    !>
    !>  \param[in]  lb      :   The input scalar of type `real` of kind \RKH, containing the lower limit of integration.<br>
    !>                          (**optional**, default = `0`)<br>
    !>  \param[in]  ub      :   The input scalar of the same type and kind as `lb`, containing the upper limit of integration.<br>
    !>                          (**optional**, default = `3.`)<br>
    !>
    !>  \interface{int5_type}
    !>  \code{.F90}
    !>
    !>      use pm_quadTest, only: int5_type
    !>      type(int5_type) :: integrand
    !>
    !>      integrand = int5_type(lb = lb, ub = ub)
    !>      print *, "description: ", integrand%desc
    !>      print *, "lower limit: ", integrand%lb
    !>      print *, "upper limit: ", integrand%ub
    !>      print *, "Example integrand value: ", integrand%get(x)
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [integrand_type](@ref pm_quadTest::integrand_type)<br>
    !>
    !>  \test
    !>  [test_pm_quadPack](@ref test_pm_quadPack)
    !>
    !>  \final{int5_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    type, extends(integrand_type)   :: int5_type
    contains
        procedure                   :: get => getInt5
    end type

    !>  \cond excluded
    interface int5_type
    module function int5_typer(lb, ub) result(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: int5_typer
#endif
        use pm_kind, only: RKG => RKH
        real(RKG), intent(in)   :: lb, ub
        type(int5_type)         :: self
    end function
    end interface
    !>  \endcond excluded

    interface
    module function getInt5(self, x) result(func)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInt5
#endif
        use pm_kind, only: RKG => RKH
        class(int5_type), intent(in)    :: self
        real(RKG)       , intent(in)    :: x
        real(RKG)                       :: func
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for generating test integrand objects of the following algebraic form.
    !>
    !>  \details
    !>  The full integrand is defined as,
    !>  \f{equation}{
    !>
    !>      f(x) = \frac{\log(x)}{1 + 100x^2} , x \in (0, +\infty)
    !>
    !>  \f}
    !>  with an integral of \f$-\pi\frac{\log(10)}{20}\f$.<br>
    !>  This integrand is inspired by the examples of John Burkardt test suite for QAGI routine of QuadPack.<br>
    !>
    !>  \interface{int6_type}
    !>  \code{.F90}
    !>
    !>      use pm_quadTest, only: int6_type
    !>      type(int6_type) :: integrand
    !>
    !>      integrand = int6_type()
    !>      print *, "description: ", integrand%desc
    !>      print *, "lower limit: ", integrand%lb
    !>      print *, "upper limit: ", integrand%ub
    !>      print *, "Example integrand value: ", integrand%get(x)
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [integrand_type](@ref pm_quadTest::integrand_type)<br>
    !>
    !>  \test
    !>  [test_pm_quadPack](@ref test_pm_quadPack)
    !>
    !>  \final{int6_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    type, extends(integrand_type)   :: int6_type
    contains
        procedure                   :: get => getInt6
    end type

    !>  \cond excluded
    interface int6_type
    module function int6_typer() result(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: int6_typer
#endif
        use pm_kind, only: RKG => RKH
        type(int6_type) :: self
    end function
    end interface
    !>  \endcond excluded

    interface
    module function getInt6(self, x) result(func)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInt6
#endif
        use pm_kind, only: RKG => RKH
        class(int6_type), intent(in)    :: self
        real(RKG)       , intent(in)    :: x
        real(RKG)                       :: func
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for generating test integrand objects of the following algebraic form.
    !>
    !>  \details
    !>  The full integrand is defined as,
    !>  \f{equation}{
    !>
    !>      f(x) = \frac{\log\big( |(1 - x^2)(1 - 2x^2)| \big) - 4\log(x)}{x^5} , x \in (\frac{1}{3}, +\infty)
    !>
    !>  \f}
    !>  with an integral of \f$52.7407483834714449977291997202299809\f$.<br>
    !>  The integrand has singularities and break points at `break = [1 / sqrt(2), 1]`.<br>
    !>
    !>  \interface{int7_type}
    !>  \code{.F90}
    !>
    !>      use pm_quadTest, only: int7_type
    !>      type(int7_type) :: integrand
    !>
    !>      integrand = int7_type()
    !>      print *, "description: ", integrand%desc
    !>      print *, "lower limit: ", integrand%lb
    !>      print *, "upper limit: ", integrand%ub
    !>      print *, "Example integrand value: ", integrand%get(x)
    !>
    !>  \endcode
    !>
    !>  \remark
    !>  This integrand is a transformation of [int5_type](@ref pm_quadTest::int5_type) such that
    !>  the integral of `int5_type%get(1/x) / x**2` over `(0, 3)` equals
    !>  the integral of `int7_type%get(x)` over `(1/3, +Inf)`.<br>
    !>
    !>  \see
    !>  [integrand_type](@ref pm_quadTest::integrand_type)<br>
    !>
    !>  \test
    !>  [test_pm_quadPack](@ref test_pm_quadPack)
    !>
    !>  \final{int7_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    type, extends(integrand_type)   :: int7_type
    contains
        procedure                   :: get => getInt7
    end type

    !>  \cond excluded
    interface int7_type
    module function int7_typer() result(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: int7_typer
#endif
        use pm_kind, only: RKG => RKH
        type(int7_type) :: self
    end function
    end interface
    !>  \endcond excluded

    interface
    module function getInt7(self, x) result(func)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInt7
#endif
        use pm_kind, only: RKG => RKH
        class(int7_type), intent(in)    :: self
        real(RKG)       , intent(in)    :: x
        real(RKG)                       :: func
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for generating test integrand objects of the following algebraic form.
    !>
    !>  \details
    !>  The full integrand is defined as,
    !>  \f{equation}{
    !>
    !>      f(x) = -\frac{\log\big( |(1 - x^2)(1 - 2x^2)| \big) - 4\log(x)}{x^5} , x \in (-\infty, -\frac{1}{3})
    !>
    !>  \f}
    !>  with an integral of \f$52.7407483834714449977291997202299809\f$.<br>
    !>  The integrand has singularities and break points at `break = [-1, -1 / sqrt(2)]`.<br>
    !>
    !>  \interface{int8_type}
    !>  \code{.F90}
    !>
    !>      use pm_quadTest, only: int8_type
    !>      type(int8_type) :: integrand
    !>
    !>      integrand = int8_type()
    !>      print *, "description: ", integrand%desc
    !>      print *, "lower limit: ", integrand%lb
    !>      print *, "upper limit: ", integrand%ub
    !>      print *, "Example integrand value: ", integrand%get(x)
    !>
    !>  \endcode
    !>
    !>  \remark
    !>  This integrand is a transformation of [int7_type](@ref pm_quadTest::int7_type) such that
    !>  the integral of `-int7_type%get(-x)` over `(1/3, +Inf)` equals
    !>  the integral of `int8_type%get(x)` over `(-Inf, -1/3)`.<br>
    !>
    !>  \see
    !>  [integrand_type](@ref pm_quadTest::integrand_type)<br>
    !>
    !>  \test
    !>  [test_pm_quadPack](@ref test_pm_quadPack)
    !>
    !>  \final{int8_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    type, extends(integrand_type)   :: int8_type
    contains
        procedure                   :: get => getInt8
    end type

    !>  \cond excluded
    interface int8_type
    module function int8_typer() result(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: int8_typer
#endif
        use pm_kind, only: RKG => RKH
        type(int8_type) :: self
    end function
    end interface
    !>  \endcond excluded

    interface
    module function getInt8(self, x) result(func)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInt8
#endif
        use pm_kind, only: RKG => RKH
        class(int8_type), intent(in)    :: self
        real(RKG)       , intent(in)    :: x
        real(RKG)                       :: func
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for generating test integrand objects of the following algebraic form.
    !>
    !>  \details
    !>  The full integrand is defined as,
    !>  \f{equation}{
    !>
    !>      f(x) =
    !>      \begin{cases}
    !>          \frac{\log\big( |(1 - x^2)(1 - 2x^2)| \big) - 4\log(x)}{x^5} &, x \in (\frac{1}{3}, +\infty) \\
    !>          \frac{1}{\pi\sqrt{-(9 + x)(10 + x)}} &, x \in (-10, -9) \\
    !>          0 &, \text{otherwise}
    !>      \end{cases}
    !>
    !>  \f}
    !>  with an integral of \f$53.7407483834714449977291997202299809\f$.<br>
    !>  The integrand has singularities and break points at `break = [-10, -9, 1 / sqrt(2), 1]`.<br>
    !>
    !>  \interface{int9_type}
    !>  \code{.F90}
    !>
    !>      use pm_quadTest, only: int9_type
    !>      type(int9_type) :: integrand
    !>
    !>      integrand = int9_type()
    !>      print *, "description: ", integrand%desc
    !>      print *, "lower limit: ", integrand%lb
    !>      print *, "upper limit: ", integrand%ub
    !>      print *, "Example integrand value: ", integrand%get(x)
    !>
    !>  \endcode
    !>
    !>  \remark
    !>  This integrand is a transformation of [int7_type](@ref pm_quadTest::int7_type)
    !>  mixed with the PDF of the Beta distribution for the shape parameters \f$(0.5, 0.5)\f$.<br>
    !>
    !>  \see
    !>  [integrand_type](@ref pm_quadTest::integrand_type)<br>
    !>
    !>  \test
    !>  [test_pm_quadPack](@ref test_pm_quadPack)
    !>
    !>  \final{int9_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    type, extends(integrand_type)   :: int9_type
    contains
        procedure                   :: get => getInt9
    end type

    !>  \cond excluded
    interface int9_type
    module function int9_typer() result(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: int9_typer
#endif
        use pm_kind, only: RKG => RKH
        type(int9_type) :: self
    end function
    end interface
    !>  \endcond excluded

    interface
    module function getInt9(self, x) result(func)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInt9
#endif
        use pm_kind, only: RKG => RKH
        class(int9_type), intent(in)    :: self
        real(RKG)       , intent(in)    :: x
        real(RKG)                       :: func
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for generating test integrand objects of the following algebraic form.
    !>
    !>  \details
    !>  The full integrand is defined as,
    !>  \f{equation}{
    !>
    !>      \large
    !>      f(x) = \left(\frac{x}{\ms{lb}} \right)^\alpha \exp\left( -\beta \left[ x - \ms{lb} \right] \right) ~, \ms{lb} \in (0, +\infty)
    !>
    !>  \f}
    !>  where \f$\beta > 0\f$ with integration range as \f$[\ms{lb}, \ms{ub}]\f$ where \f$\ms{lb} < \ms{ub} < +\infty\f$.<br>
    !>  The integrand has a singularity at \f$x = 0\f$ with \f$\alpha < 0\f$, but the \f$\ms{lb}\f$ range does not allow singularity to enter the integrand.<br>
    !>  When \f$\alpha > -1, \ms{ub} = +\infty\f$, this integral can be computed via [regularized upper incomplete Gamma function](@ref pm_mathGamma) \f$Q(\cdot)\f$:<br>
    !>  \f{equation}{
    !>
    !>      \large
    !>      f(x) 
    !>      = \frac{ \exp(\beta \ms{lb}) }{ \ms{lb}^\alpha ~ \beta^{\alpha + 1} } ~ \Gamma(\alpha + 1) ~ Q(\alpha + 1, \beta \ms{lb})
    !>      - \frac{ \exp(\beta \ms{ub}) }{ \ms{ub}^\alpha ~ \beta^{\alpha + 1} } ~ \Gamma(\alpha + 1) ~ Q(\alpha + 1, \beta \ms{ub})
    !>      ~.
    !>
    !>  \f}
    !>  Otherwise, the integrand must be computed numerically, in which case, the `integrand` component of object (representing the truth) is set to `NaN`.<br>
    !>
    !>  \param[in]  lb      :   The input scalar of type `real` of kind \RKH.<br>
    !>                          (**optional**, default = `1.`)<br>
    !>  \param[in]  ub      :   The input scalar of the same type and kind as `a`.<br>
    !>                          (**optional**, default = [getInfPos(self%ub)](@ref pm_except::getInfPos)<br>
    !>  \param[in]  alpha   :   The input scalar of type `integer` of default kind \IK, standing for Lower Factor, such that `lb = lf * pi` is the lower bound of integration.<br>
    !>                          (**optional**, default = `+1.`)<br>
    !>  \param[in]  beta    :   The input scalar of type `integer` of default kind \IK, standing for Upper Factor, such that `ub = uf * pi` is the upper bound of integration.<br>
    !>                          (**optional**, default = `+1.`)<br>
    !>
    !>  \interface{intGamUpp_type}
    !>  \code{.F90}
    !>
    !>      use pm_quadTest, only: intGamUpp_type
    !>      type(intGamUpp_type) :: integrand
    !>
    !>      integrand = intGamUpp_type()
    !>      print *, "description: ", integrand%desc
    !>      print *, "lower limit: ", integrand%alpha
    !>      print *, "lower limit: ", integrand%beta
    !>      print *, "lower limit: ", integrand%lb
    !>      print *, "upper limit: ", integrand%ub
    !>      print *, "Example integrand value: ", integrand%get(x)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < lb` must hold for the corresponding input arguments.<br>
    !>  The condition `lb < ub` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < beta` must hold for the corresponding input arguments.<br>
    !>
    !>  \see
    !>  [integrand_type](@ref pm_quadTest::integrand_type)<br>
    !>
    !>  \test
    !>  [test_pm_quadPack](@ref test_pm_quadPack)
    !>
    !>  \final{intGamUpp_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    type, extends(integrand_type)   :: intGamUpp_type
        real(RKH)                   :: alpha, beta
        real(RKH)                   :: normfac
    contains
        procedure                   :: get => getIntGamUpp
    end type

    !>  \cond excluded
    interface intGamUpp_type
    module function intGamUpp_typer(lb, ub, alpha, beta) result(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: intGamUpp_typer
#endif
        use pm_kind, only: RKG => RKH
        real(RKG), intent(in), optional :: lb, ub, alpha, beta
        type(intGamUpp_type) :: self
    end function
    end interface
    !>  \endcond excluded

    interface
    module function getIntGamUpp(self, x) result(func)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getIntGamUpp
#endif
        use pm_kind, only: RKG => RKH
        class(intGamUpp_type), intent(in) :: self
        real(RKG), intent(in) :: x
        real(RKG) :: func
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for generating test integrand objects of the trigonometric form as described below.
    !>
    !>  \details
    !>  The full integrand is defined over a finite interval as,<br>
    !>  \f{equation}{
    !>
    !>      f(x) = \cos(a\sin(bx)) ~,~ x \in (-\infty < \ms{lb} = \mathrm{lf} * \pi, \ms{ub} = \mathrm{uf} * \pi < +\infty)
    !>
    !>  \f}
    !>  where the factors \f$a\f$ and \f$b\f$ are any finite real numbers and \f$(\mathrm{lf}, \mathrm{uf})\f$ are whole numbers (integer-valued).<br>
    !>  The definite integral of the integrand is,<br>
    !>  \f{equation}{
    !>
    !>      \int_{\ms{lb}}^{\ms{ub}} f(x) dx = (\ms{ub} - \ms{lb}) J_0(a) ~,
    !>
    !>  \f}
    !>  where \f$J_0\f$ is the Modified Bessel function of the zeroth kind.<br>
    !>
    !>  \param[in]  lf      :   The input scalar of type `integer` of default kind \IK, standing for Lower Factor, such that `lb = lf * pi` is the lower bound of integration.<br>
    !>                          (**optional**, default = `-1.`)<br>
    !>  \param[in]  uf      :   The input scalar of type `integer` of default kind \IK, standing for Upper Factor, such that `ub = uf * pi` is the upper bound of integration.<br>
    !>                          (**optional**, default = `+1.`)<br>
    !>  \param[in]  a       :   The input scalar of type `real` of kind \RKH.<br>
    !>                          (**optional**, default = `10.`)<br>
    !>  \param[in]  b       :   The input scalar of the same type and kind as `a`.<br>
    !>                          (**optional**, default = `+1.`)<br>
    !>
    !>  \interface{intSinCos_type}
    !>  \code{.F90}
    !>
    !>      use pm_quadTest, only: intSinCos_type
    !>      type(intSinCos_type) :: integrand
    !>
    !>      integrand = intSinCos_type(lb = lb, ub = ub, a = a, b = b)
    !>      print *, "description: ", integrand%desc
    !>      print *, "lower limit: ", integrand%lb
    !>      print *, "upper limit: ", integrand%ub
    !>      print *, "Example integrand value: ", integrand%get(x)
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [integrand_type](@ref pm_quadTest::integrand_type)<br>
    !>
    !>  \test
    !>  [test_pm_quadPack](@ref test_pm_quadPack)
    !>
    !>  \final{intSinCos_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    type, extends(integrand_type)   :: intSinCos_type
        integer(IK)                 :: lf, uf
        real(RKH)                   :: a, b
    contains
        procedure                   :: get => getIntSinCos
    end type

    !>  \cond excluded
    interface intSinCos_type
    module function intSinCos_typer(lf, uf, a, b) result(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: intSinCos_typer
#endif
        use pm_kind, only: RKG => RKH
        integer(IK) , intent(in), optional  :: lf, uf
        real(RKG)   , intent(in), optional  :: a, b
        type(intSinCos_type)                :: self
    end function
    end interface
    !>  \endcond excluded

    interface
    module function getIntSinCos(self, x) result(func)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getIntSinCos
#endif
        use pm_kind, only: RKG => RKH
        class(intSinCos_type)  , intent(in) :: self
        real(RKG)           , intent(in)    :: x
        real(RKG)                           :: func
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for generating test integrand objects of the Probability Density Function of the Normal distribution.
    !>
    !>  \details
    !>  The full integrand is defined as,
    !>  \f{equation}{
    !>
    !>  \pi(x | \mu, \sigma) = \frac{1}{\sigma\sqrt{2\pi}}\exp\bigg( -\frac{\big(x - \mu\big)^2}{2\sigma^2} \bigg) ~,~ x \in (-\infty, +\infty)
    !>
    !>  \f}
    !>
    !>  \param[in]  lb      :   The input scalar of type `real` of kind \RKH, containing the lower limit of integration.<br>
    !>                          (**optional**, default = [getInfNeg(real(0,kind(lb))](@ref pm_except::getInfNeg))<br>
    !>  \param[in]  ub      :   The input scalar of the same type and kind as `lb`, containing the upper limit of integration.<br>
    !>                          (**optional**, default = [getInfPos(real(0,kind(ub))](@ref pm_except::getInfPos))<br>
    !>  \param[in]  mu      :   The input scalar of the same type and kind as `lb`, representing the location parameter of the Normal distribution.<br>
    !>                          (**optional**, default = `0`)<br>
    !>  \param[in]  sigma   :   The input scalar of the same type and kind as `lb`, representing the scale parameter of the Normal distribution.<br>
    !>                          (**optional**, default = `1.`)<br>
    !>
    !>  \interface{intNormPDF_type}
    !>  \code{.F90}
    !>
    !>      use pm_quadTest, only: intNormPDF_type
    !>      type(intNormPDF_type) :: integrand
    !>
    !>      integrand = intNormPDF_type(lb = lb, ub = ub, mu = mu, sigma = sigma)
    !>      print *, "description: ", integrand%desc
    !>      print *, "lower limit: ", integrand%lb
    !>      print *, "upper limit: ", integrand%ub
    !>      print *, "Example integrand value: ", integrand%get(x)
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [integrand_type](@ref pm_quadTest::integrand_type)<br>
    !>
    !>  \test
    !>  [test_pm_quadPack](@ref test_pm_quadPack)
    !>
    !>  \final{intNormPDF_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    type, extends(integrand_type)   :: intNormPDF_type
        real(RKH)                   :: mu               !<  \public The location parameter of the Normal distribution.
        real(RKH)                   :: sigma            !<  \public The scale parameter (standard deviation) of the Normal distribution.
        real(RKH)                   :: invSigma         !<  \public The inverse scale parameter (standard deviation) of the Normal distribution.
        real(RKH)                   :: logInvSigma      !<  \public The natural logarithm of the inverse scale parameter (standard deviation) of the Normal distribution.
    contains
        procedure                   :: get => getIntNormPDF
    end type

    !>  \cond excluded
    interface intNormPDF_type
    module function intNormPDF_typer(lb, ub, mu, sigma) result(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: intNormPDF_typer
#endif
        use pm_kind, only: RKG => RKH
        real(RKG), intent(in), optional :: lb, ub, mu, sigma
        type(intNormPDF_type)           :: self
    end function
    end interface
    !>  \endcond excluded

    interface
    module function getIntNormPDF(self, x) result(func)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getIntNormPDF
#endif
        use pm_kind, only: RKG => RKH
        class(intNormPDF_type)  , intent(in)    :: self
        real(RKG)               , intent(in)    :: x
        real(RKG)                               :: func
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for generating test integrand objects of the Probability Density Function of the Lognormal distribution.
    !>
    !>  \details
    !>  The full integrand is defined as,
    !>  \f{equation}{
    !>
    !>  \pi(x | \mu, \sigma) = \frac{1}{x\sigma\sqrt{2\pi}}\exp\bigg( -\frac{\big(\log(x) - \mu\big)^2}{2\sigma^2} \bigg) ~,~ x \in (-\infty, +\infty)
    !>
    !>  \f}
    !>  with an integral value of `1`.
    !>
    !>  \param[in]  lb      :   The input scalar of type `real` of kind \RKH, containing the lower limit of integration.<br>
    !>                          (**optional**, default = `0`)<br>
    !>  \param[in]  ub      :   The input scalar of the same type and kind as `lb`, containing the upper limit of integration.<br>
    !>                          (**optional**, default = [getInfPos(real(0,kind(ub))](@ref pm_except::getInfPos))<br>
    !>  \param[in]  mu      :   The input scalar of the same type and kind as `lb`, representing the location parameter of the Lognormal distribution.<br>
    !>                          (**optional**, default = `0`)<br>
    !>  \param[in]  sigma   :   The input scalar of the same type and kind as `lb`, representing the scale parameter of the Lognormal distribution.<br>
    !>                          (**optional**, default = `1.`)<br>
    !>
    !>  \interface{intLogNormPDF_type}
    !>  \code{.F90}
    !>
    !>      use pm_quadTest, only: intLogNormPDF_type
    !>      type(intLogNormPDF_type) :: integrand
    !>
    !>      integrand = intLogNormPDF_type(lb = lb, ub = ub, mu = mu, sigma = sigma)
    !>      print *, "description: ", integrand%desc
    !>      print *, "lower limit: ", integrand%lb
    !>      print *, "upper limit: ", integrand%ub
    !>      print *, "Example integrand value: ", integrand%get(x)
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [integrand_type](@ref pm_quadTest::integrand_type)<br>
    !>
    !>  \test
    !>  [test_pm_quadPack](@ref test_pm_quadPack)
    !>
    !>  \final{intLogNormPDF_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    type, extends(integrand_type)   :: intLogNormPDF_type
        real(RKH)                   :: mu               !<  \public The location parameter of the Normal distribution.
        real(RKH)                   :: sigma            !<  \public The scale parameter (standard deviation) of the Lognormal distribution.
        real(RKH)                   :: invSigma         !<  \public The inverse scale parameter (standard deviation) of the Lognormal distribution.
        real(RKH)                   :: logInvSigma      !<  \public The natural logarithm of the inverse scale parameter (standard deviation) of the Lognormal distribution.
    contains
        procedure                   :: get => getIntLogNormPDF
    end type

    !>  \cond excluded
    interface intLogNormPDF_type
    module function intLogNormPDF_typer(lb, ub, mu, sigma) result(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: intLogNormPDF_typer
#endif
        use pm_kind, only: RKG => RKH
        real(RKG), intent(in), optional :: lb, ub, mu, sigma
        type(intLogNormPDF_type)        :: self
    end function
    end interface
    !>  \endcond excluded

    interface
    module function getIntLogNormPDF(self, x) result(func)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getIntLogNormPDF
#endif
        use pm_kind, only: RKG => RKH
        class(intLogNormPDF_type)   , intent(in)    :: self
        real(RKG)                   , intent(in)    :: x
        real(RKG)                                   :: func
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for generating test integrand objects of the Probability Density Function of the GenExpGamma distribution.
    !>
    !>  \details
    !>  The full integrand is the PDF of the GenExpGamma distribution as defined in the documentation of [pm_distGenExpGamma](@ref pm_distGenExpGamma),
    !>  with an integral value of `1` over its full support \f$x \in (-\infty, +\infty)\f$.
    !>
    !>  \param[in]  lb          :   The input scalar of type `real` of kind \RKH, containing the lower limit of integration.<br>
    !>                              (**optional**, default = [getInfNeg(real(0,kind(lb))](@ref pm_except::getInfNeg))<br>
    !>  \param[in]  ub          :   The input scalar of the same type and kind as `lb`, containing the upper limit of integration.<br>
    !>                              (**optional**, default = [getInfPos(real(0,kind(ub))](@ref pm_except::getInfPos))<br>
    !>  \param[in]  kappa       :   The input scalar of the same type and kind as `lb`, representing the shape parameter of the GenExpGamma distribution.<br>
    !>                              (**optional**, default = `1.`)<br>
    !>  \param[in]  invOmega    :   The input scalar of the same type and kind as `lb`, representing the inverse of the scale parameter of the GenExpGamma distribution.<br>
    !>                              (**optional**, default = `1.`)<br>
    !>  \param[in]  logSigma    :   The input scalar of the same type and kind as `lb`, representing the location parameter of the GenExpGamma distribution.<br>
    !>                              (**optional**, default = `0`)<br>
    !>
    !>  \interface{intGenExpGammaPDF_type}
    !>  \code{.F90}
    !>
    !>      use pm_quadTest, only: intGenExpGammaPDF_type
    !>      type(intGenExpGammaPDF_type) :: integrand
    !>
    !>      integrand = intGenExpGammaPDF_type(lb = lb, ub = ub, kappa = kappa, invOmega = invOmega, logSigma = logSigma)
    !>      print *, "description: ", integrand%desc
    !>      print *, "lower limit: ", integrand%lb
    !>      print *, "upper limit: ", integrand%ub
    !>      print *, "Example integrand value: ", integrand%get(x)
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [integrand_type](@ref pm_quadTest::integrand_type)<br>
    !>
    !>  \test
    !>  [test_pm_quadPack](@ref test_pm_quadPack)
    !>
    !>  \final{intGenExpGammaPDF_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    type, extends(integrand_type)   :: intGenExpGammaPDF_type
        real(RKH)                   :: kappa        !<  \public The shape parameter of the GenExpGamma distribution.
        real(RKH)                   :: invOmega     !<  \public The inverse of the scale parameter of the GenExpGamma distribution.
        real(RKH)                   :: logSigma     !<  \public The location parameter of the GenExpGamma distribution.
        real(RKH)                   :: logPDFNF     !<  \public The natural logarithm of the normalization factor of the GenExpGamma distribution.
    contains
        procedure                   :: get => getIntGenExpGammaPDF
    end type

    !>  \cond excluded
    interface intGenExpGammaPDF_type
    module function intGenExpGammaPDF_typer(lb, ub, kappa, invOmega, logSigma) result(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: intGenExpGammaPDF_typer
#endif
        use pm_kind, only: RKG => RKH
        real(RKG)   , intent(in), optional  :: lb, ub, kappa, invOmega, logSigma
        type(intGenExpGammaPDF_type)        :: self
    end function
    end interface
    !>  \endcond excluded

    interface
    module function getIntGenExpGammaPDF(self, x) result(func)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getIntGenExpGammaPDF
#endif
        use pm_kind, only: RKG => RKH
        class(intGenExpGammaPDF_type)   , intent(in)    :: self
        real(RKG)                       , intent(in)    :: x
        real(RKG)                                       :: func
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for generating test integrand objects of the sum of five Probability Density Functions of the Gamma distribution.
    !>
    !>  \details
    !>  The full integrand is defined as,
    !>  \f{eqnarray}{
    !>
    !>      f(x)    &=& \pi_\mathcal{G}(x + 9; 0.7, 1) \\ \nonumber
    !>              &+& \pi_\mathcal{G}(x + 5; 0.7, 1) \\ \nonumber
    !>              &+& \pi_\mathcal{G}(x - 5; 0.7, 1) \\ \nonumber
    !>              &+& \pi_\mathcal{G}(2 - x; 0.7, 1) \\ \nonumber
    !>              &+& \pi_\mathcal{G}(7 - x; 0.7, 1) \\ \nonumber
    !>              &,& x \in (-\infty, +\infty)
    !>
    !>  \f}
    !>  where \f$\pi_\mathcal{G}\f$ represents the PDF of the Gamma distribution computed via [getGammaLogPDF](@ref pm_distGamma::getGammaLogPDF).<br>
    !>  The integrand has five singularities and break points at \f$\ms{break} = [-9, -5, 2, 5, 7]\f$.<br>
    !>  By definition, the integral of the integrand over the entire fully-infinite integration range is `5.`.<br>
    !>
    !>  \interface{intPentaGammaInf_type}
    !>  \code{.F90}
    !>
    !>      use pm_quadTest, only: intPentaGammaInf_type
    !>      type(intPentaGammaInf_type) :: integrand
    !>
    !>      integrand = intPentaGammaInf_type()
    !>      print *, "description: ", integrand%desc
    !>      print *, "lower limit: ", integrand%lb
    !>      print *, "upper limit: ", integrand%ub
    !>      print *, "Example integrand value: ", integrand%get(x)
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [integrand_type](@ref pm_quadTest::integrand_type)<br>
    !>
    !>  \test
    !>  [test_pm_quadPack](@ref test_pm_quadPack)
    !>
    !>  \final{intPentaGammaInf_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    type, extends(integrand_type)   :: intPentaGammaInf_type
    contains
        procedure                   :: get => getIntPentaGammaInf
    end type

    !>  \cond excluded
    interface intPentaGammaInf_type
    module function intPentaGammaInf_typer() result(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: intPentaGammaInf_typer
#endif
        use pm_kind, only: RKG => RKH
        type(intPentaGammaInf_type) :: self
    end function
    end interface
    !>  \endcond excluded

    interface
    module function getIntPentaGammaInf(self, x) result(func)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getIntPentaGammaInf
#endif
        use pm_kind, only: RKG => RKH
        class(intPentaGammaInf_type), intent(in)    :: self
        real(RKG)                   , intent(in)    :: x
        real(RKG)                                   :: func
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for generating test integrand objects of algebraic form as described below.
    !>
    !>  \details
    !>  The full integrand is defined over a finite interval as,<br>
    !>  \f{equation}{
    !>
    !>      f(x) = \frac{1}{(1 + x)\sqrt{x}} ~,~ x \in (0 < \ms{lb}, \ms{ub}) ~.
    !>
    !>  \f}
    !>  The integrand has the precise value,<br>
    !>  \f{equation}{
    !>
    !>      \int_{\ms{lb}}^{\ms{ub}} f(x) dx = 2 ( \mathrm{atan}(\sqrt{\ms{ub}}) - \mathrm{atan}(\sqrt{\ms{lb}}) ) ~.
    !>
    !>  \f}
    !>  This integrand is an extension of the example discussed in Doncker et al (1976), Automatic Computation of Integrals with Singular integrand.<br>
    !>
    !>  \param[in]  lb      :   The input negative scalar of type `real` of kind \RKH, containing the lower limit of integration.<br>
    !>                          (**optional**, default = `0._RK`)<br>
    !>  \param[in]  ub      :   The input positive scalar of the same type and kind as `lb`, containing the upper limit of integration.<br>
    !>                          (**optional**, default = [getInfPos(0._RKG)](@ref pm_except::getInfPos))<br>
    !>
    !>  \interface{intDoncker1_type}
    !>  \code{.F90}
    !>
    !>      use pm_quadTest, only: intDoncker1_type
    !>      type(intSinCos_type) :: integrand
    !>
    !>      integrand = intDoncker1_type(lb = lb, ub = ub)
    !>      print *, "description: ", integrand%desc
    !>      print *, "lower limit: ", integrand%lb
    !>      print *, "upper limit: ", integrand%ub
    !>      print *, "Example integrand value: ", integrand%get(x)
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [integrand_type](@ref pm_quadTest::integrand_type)<br>
    !>
    !>  \test
    !>  [test_pm_quadPack](@ref test_pm_quadPack)
    !>
    !>  \final{intDoncker1_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    type, extends(integrand_type)   :: intDoncker1_type
    contains
        procedure                   :: get => getIntDoncker1
    end type

    !>  \cond excluded
    interface intDoncker1_type
    module function intDoncker1_typer(lb, ub) result(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: intDoncker1_typer
#endif
        use pm_kind, only: RKG => RKH
        real(RKG), intent(in), optional :: lb, ub
        type(intDoncker1_type)          :: self
    end function
    end interface
    !>  \endcond excluded

    interface
    module function getIntDoncker1(self, x) result(func)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getIntDoncker1
#endif
        use pm_kind, only: RKG => RKH
        class(intDoncker1_type) , intent(in)    :: self
        real(RKG)               , intent(in)    :: x
        real(RKG)                               :: func
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for generating test integrand objects of algebraic form as described below.
    !>
    !>  \details
    !>  The full integrand is defined over a finite interval as,<br>
    !>  \f{equation}{
    !>
    !>      f(x) = \frac{\exp(x)}{\sqrt{-x}} ~,~ x \in (\ms{lb}, \ms{ub} < 0.) ~.
    !>
    !>  \f}
    !>  The integrand has the precise value,<br>
    !>  \f{equation}{
    !>
    !>      \int_{\ms{lb} = -\infty}^{\ms{ub}} f(x) dx = \sqrt{\pi} \big(\mathrm{erf}(\sqrt{-\ms{lb}} - \mathrm{erf}(\sqrt{-\ms{ub}}) \big) ~.
    !>
    !>  \f}
    !>  This integrand is an extension of the example discussed in Doncker et al (1976), Automatic Computation of Integrals with Singular integrand.<br>
    !>
    !>  \param[in]  lb      :   The input negative scalar of type `real` of kind \RKH, containing the lower limit of integration.<br>
    !>                          (**optional**, default = [getInfNeg(0._RKG)](@ref pm_except::getInfNeg))<br>
    !>  \param[in]  ub      :   The input positive scalar of the same type and kind as `lb`, containing the upper limit of integration.<br>
    !>                          (**optional**, default = `0._RK`)<br>
    !>
    !>  \interface{intDoncker2_type}
    !>  \code{.F90}
    !>
    !>      use pm_quadTest, only: intDoncker2_type
    !>      type(intSinCos_type) :: integrand
    !>
    !>      integrand = intDoncker2_type(lb = lb, ub = ub)
    !>      print *, "description: ", integrand%desc
    !>      print *, "lower limit: ", integrand%lb
    !>      print *, "upper limit: ", integrand%ub
    !>      print *, "Example integrand value: ", integrand%get(x)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The conditions `lb < ub < 0._RK` must for the relevant input arguments.<br>
    !>  \vericons
    !>
    !>  \see
    !>  [integrand_type](@ref pm_quadTest::integrand_type)<br>
    !>
    !>  \test
    !>  [test_pm_quadPack](@ref test_pm_quadPack)
    !>
    !>  \final{intDoncker2_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    type, extends(integrand_type)   :: intDoncker2_type
    contains
        procedure                   :: get => getIntDoncker2
    end type

    !>  \cond excluded
    interface intDoncker2_type
    module function intDoncker2_typer(lb, ub) result(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: intDoncker2_typer
#endif
        use pm_kind, only: RKG => RKH
        real(RKG), intent(in), optional :: lb, ub
        type(intDoncker2_type)          :: self
    end function
    end interface
    !>  \endcond excluded

    interface
    module function getIntDoncker2(self, x) result(func)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getIntDoncker2
#endif
        use pm_kind, only: RKG => RKH
        class(intDoncker2_type) , intent(in)    :: self
        real(RKG)               , intent(in)    :: x
        real(RKG)                               :: func
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for generating test integrand objects of the algebraic form as described below, whose Cauchy Principal Value is to be computed.
    !>
    !>  \details
    !>  The full integrand is defined as,<br>
    !>  \f{equation}{
    !>      f(x) = \frac{1}{(x - c)} ~,~ x \in (\ms{lb}, \ms{ub}) ~,~ \ms{lb} < 0 < \ms{ub}
    !>  \f}
    !>  where the integration bounds are finite values.<br>
    !>  The Cauchy Principal value of the integrand is \f$\log\bigg(\frac{\ms{ub} - \ms{cs}}{\ms{lb} - \ms{cs}}\bigg)\f$.<br>
    !>
    !>  \param[in]  lb      :   The input scalar of type `real` of kind \RKH, containing the lower limit of integration.<br>
    !>                          (**optional**, default = `-2.`)<br>
    !>  \param[in]  ub      :   The input scalar of the same type and kind as `lb`, containing the upper limit of integration.<br>
    !>                          (**optional**, default = `+3.`)<br>
    !>  \param[in]  cs      :   The input scalar of the same type and kind as `lb`, containing the Cauchy singularity of the integrand.<br>
    !>                          (**optional**, default = `+1.`)<br>
    !>
    !>  \interface{intCauchy1_type}
    !>  \code{.F90}
    !>
    !>      use pm_quadTest, only: intCauchy1_type
    !>      type(intCauchy1_type) :: integrand
    !>
    !>      integrand = intCauchy1_type(lb = lb, ub = ub, cs = cs)
    !>      print *, "description: ", integrand%desc
    !>      print *, "lower limit: ", integrand%lb
    !>      print *, "upper limit: ", integrand%ub
    !>      print *, "singularity: ", integrand%cs
    !>      print *, "Example integrand value: ", integrand%get(x)
    !>      print *, "Example integrand value without the Cauchy weight: ", integrand%getWeighted(x)
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [integrand_type](@ref pm_quadTest::integrand_type)<br>
    !>
    !>  \test
    !>  [test_pm_quadPack](@ref test_pm_quadPack)
    !>
    !>  \final{intCauchy1_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    type, extends(integrand_type)   :: intCauchy1_type
    contains
        procedure                   :: get => getIntCauchy1
    end type

    !>  \cond excluded
    interface intCauchy1_type
    module function intCauchy1_typer(lb, ub, cs) result(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: intCauchy1_typer
#endif
        use pm_kind, only: RKG => RKH
        real(RKG)   , intent(in), optional  :: lb, ub, cs
        type(intCauchy1_type)               :: self
    end function
    end interface
    !>  \endcond excluded

    interface
    module function getIntCauchy1(self, x) result(func)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getIntCauchy1
#endif
        use pm_kind, only: RKG => RKH
        class(intCauchy1_type)  , intent(in)    :: self
        real(RKG)               , intent(in)    :: x
        real(RKG)                               :: func
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type for generating test integrand objects of the algebraic form as described below, whose Cauchy Principal Value is to be computed.
    !>
    !>  \details
    !>  The full integrand is defined as,<br>
    !>  \f{equation}{
    !>
    !>      f(x) = \frac{1}{(x - \ms{cs1})(x - \ms{cs2})}
    !>      ~,~
    !>      x \in (-\infty \leq \ms{lb} < \min(\ms{cs1},\ms{cs2}), \min(\ms{cs1},\ms{cs2}) < \ms{ub} < \max(\ms{cs1},\ms{cs2}))
    !>      ~\vee~
    !>      x \in (\min(\ms{cs1},\ms{cs2}) < \ms{lb} < \max(\ms{cs1},\ms{cs2}), \max(\ms{cs1},\ms{cs2}) < \ms{ub} \leq +\infty)
    !>
    !>  \f}
    !>  Depending on the choice of integration range, the integrand has either \f$\ms{cs1}\f$ or \f$\ms{cs2}\f$ as its Cauchy singularity (but not both).<br>
    !>  The Cauchy Principal value of the integrand is,<br>
    !>  \f{equation}{
    !>
    !>      \bigg[ \frac{\log(x - \ms{cs1}) - \log(x - \ms{cs2})}{\ms{cs1} - \ms{cs2}} \bigg]_{\ms{lb}}^{\ms{ub}} ~.
    !>
    !>  \f}
    !>
    !>  \param[in]  lb      :   The input scalar of type `real` of kind \RKH, containing the lower limit of integration.<br>
    !>                          (**optional**, default = `-2.`)<br>
    !>  \param[in]  ub      :   The input scalar of the same type and kind as `lb`, containing the upper limit of integration.<br>
    !>                          (**optional**, default = `+2.`)<br>
    !>  \param[in]  cs1     :   The input scalar of the same type and kind as `lb`, containing the first pole (Cauchy singularity) of the integrand.<br>
    !>                          Note that `cs1 < cs2` must hold.<br>
    !>                          (**optional**, default = `-2.`)<br>
    !>  \param[in]  cs2     :   The input scalar of the same type and kind as `lb`, containing the second pole (Cauchy singularity) of the integrand.<br>
    !>                          Note that `cs1 < cs2` must hold.<br>
    !>                          (**optional**, default = `+3.`)<br>
    !>
    !>  \interface{intCauchy2_type}
    !>  \code{.F90}
    !>
    !>      use pm_quadTest, only: intCauchy2_type
    !>      type(intCauchy2_type) :: integrand
    !>
    !>      integrand = intCauchy2_type(lb = lb, ub = ub, cs1 = cs1, cs2 = cs2)
    !>      print *, "description: ", integrand%desc
    !>      print *, "lower limit: ", integrand%lb
    !>      print *, "upper limit: ", integrand%ub
    !>      print *, "singularity: ", integrand%cs
    !>      print *, "Example integrand value: ", integrand%get(x)
    !>      print *, "Example integrand value without the Cauchy weight: ", integrand%getWeighted(x)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The two Cauchy singularities of the integrand must not be simultaneously present in the in the integration range.<br>
    !>  \vericon
    !>
    !>  \see
    !>  [integrand_type](@ref pm_quadTest::integrand_type)<br>
    !>
    !>  \test
    !>  [test_pm_quadPack](@ref test_pm_quadPack)
    !>
    !>  \final{intCauchy2_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    type, extends(integrand_type)   :: intCauchy2_type
        real(RKH)   , private       :: csnot, Pole(2)
    contains
        procedure                   :: get => getIntCauchy2
    end type

    !>  \cond excluded
    interface intCauchy2_type
    module function intCauchy2_typer(lb, ub, cs1, cs2) result(self)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: intCauchy2_typer
#endif
        use pm_kind, only: RKG => RKH
        real(RKG)   , intent(in), optional  :: lb, ub, cs1, cs2
        type(intCauchy2_type)               :: self
    end function
    end interface
    !>  \endcond excluded

    interface
    module function getIntCauchy2(self, x) result(func)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getIntCauchy2
#endif
        use pm_kind, only: RKG => RKH
        class(intCauchy2_type)  , intent(in)    :: self
        real(RKG)               , intent(in)    :: x
        real(RKG)                               :: func
    end function
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_quadTest