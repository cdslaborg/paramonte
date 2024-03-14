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
!>  This module contains the procedures and interfaces for evaluating the relative or absolute proximity of two numeric values.<br>
!>
!>  \details
!>  Floating point values contain limited precision.<br>
!>  As such, an accurate equality comparison of `real` or `complex` values is frequently impossible or insufficient.<br>
!>  Instead, one frequently needs to determine whether a computed value is **close enough** to an expected value, without requiring them to be exactly equal.<br>
!>
!>  The generic interfaces of this module implement different methods of determining whether two numeric values is sufficiently close to each other,
!>  \code{.F90}
!>      isClose(x, y, reltol = epsilon(x), abstol = tiny(x))
!>      isClose(x, y, method, reltol = epsilon(x), abstol = tiny(x))
!>  \endcode
!>  where,<br>
!>  <ul>
!>      <li>    `x` and `y` are the two values to be tested for relative closeness,<br>
!>      <li>    `method` is a scalar parameter that determines the method of computing the relative error,<br>
!>      <li>    `reltol` is the relative tolerance for the default or the specified choice of `method`,<br>
!>      <li>    `abstol` is a minimum absolute tolerance level, useful for comparisons near zero.<br>
!>  </ul>
!>
!>  \note
!>  <ul>
!>      <li>    The IEEE 754 special values of `NaN`, `inf`, and `-inf` are handled according to IEEE rules.<br>
!>      <li>    Specifically, `NaN` is **not considered close** to any other value, including `NaN`.<br>
!>      <li>    Two `+inf` values are considered close to each other at all times.<br>
!>      <li>    Two `-inf` values are considered close to each other at all times.<br>
!>      <li>    Although the generic interfaces of this module only accept `real` and `complex` values,
!>              values `integer` type can be also compared with each other by first converting them to `real` values.<br>
!>      <li>    For `complex` values, the relative and absolute tolerances are specified as `real` and the absolute value of the `complex` values will be used for scaling and comparison.<br>
!>      <li>    <b>If `complex` tolerances are desired</b>, simply pass the absolute values (`abs(reltol)` and `abs(abstol)`) of the `complex` tolerances.<br>
!>  </ul>
!>
!>  <b>Behavior near zero</b><br>
!>  Relative comparison is problematic if either value is zero.<br>
!>  Mathematically speaking, no value is small **relative** to zero.<br>
!>  Computationally speaking, if either value is zero, the difference is the absolute value of the other non-zero value,
!>  and the computed absolute tolerance will be `reltol` times that value.<br>
!>  When `reltol < 1`, the difference will never be less than the tolerance.<br>
!>  Nevertheless, there are many instances where there is a need to know if a computed value is *close to zero*.<br>
!>  This calls for *an absolute tolerance test*.<br>
!>
!>  <b>Behavior for small values straddling zero</b><br>
!>  There is a similar issue if the two values to be compared straddle zero.<br>
!>  If `x` is approximately equal to `-y`, then `x` and `y` will never be computed as *close*.<br>
!>  To handle such cases, the optional parameter, `abstol` can be used to set a minimum tolerance used in the case of very small or zero computed relative tolerance.<br>
!>  That is, <b>the values will be always be considered close if the difference between them is less than `abstol`</b>.<br>
!>  If the input `reltol` argument is set `0.`, then only the absolute tolerance will effect the result.<br>
!>  Although trivial, this allows the generic interfaces to be also used for purely absolute tolerance checks.<br>
!>
!>  <b>Relative vs. Absolute Difference</b><br>
!>  There are two ways to define the proximity of two numeric values:<br>
!>  <ol>
!>      <li>    <b>Absolute difference</b>: `abs(x - y)`,<br>
!>      <li>    <b>Relative difference</b>: `abs(x - y) / sigma`, where `sigma` is called the **scale factor**.<br>
!>              Usually, the scale factor is some function of the values under consideration, for instance:<br>
!>              <ol>
!>                  <li>    <b>Reference Scaling (asymmetric)</b>   : The absolute value of one of the input values.<br>
!>                  <li>    <b>Weak Scaling (symmetric)</b>         : The maximum absolute value of the two.<br>
!>                  <li>    <b>Strong Scaling (symmetric)</b>       : The minimum absolute value of the two.<br>
!>                  <li>    <b>mean Scaling (symmetric)</b>         : The absolute value of the arithmetic mean of the two.<br>
!>              </ol>
!>              These leads to the following possibilities for determining if two values, `x` and `y`, are close to each other.<br>
!>              <ol>
!>                  <li>    <b>Reference Scaling</b>    : `abs(x - y) <= reltol * abs(x)`.<br>
!>                  <li>    <b>Weak Scaling</b>         : `abs(x - y) <= reltol * max(abs(x), abs(y))`.<br>
!>                  <li>    <b>Strong Scaling</b>       : `abs(x - y) <= reltol * min(abs(x), abs(y))`.<br>
!>                  <li>    <b>mean Scaling</b>         : `abs(x - y) <= reltol * (x + y) / 2`.<br>
!>              </ol>
!>              The naming convention for the *weak* vs. *strong* scaling stems from the following facts:<br>
!>              <ul>
!>                  <li>    The *weak scaling* test can be also written as: `(abs(x - y) <= abs(reltol * x)) .or. (abs(x - y) <= abs(reltol * y))`.<br>
!>                  <li>    The *strong scaling* test can be also written as: `(abs(x - y) <= abs(reltol * x)) .and. (abs(x - y) <= abs(reltol * y))`.<br>
!>              </ul>
!>              Each of the above proximity testing methods can lead to slightly different results.<br>
!>              However, if `reltol` is small compared to the two input values, the differences between the methods are negligible, frequently less than available floating point precision.<br>
!>              <b>How significantly different are the relative-tolerance testing methods?</b><br>
!>              When selecting a method to determine closeness, one might want to know how much of a difference it could make to use one test or the other.<br>
!>              How many values are there (or what range of values) that will pass one test, but not the other?<br>
!>              The largest difference is between the weak and strong scaling methods above where the tolerance is scaled by either the larger or smaller of the values.<br>
!>              Suppose `delta` is the difference between the allowable absolute tolerance defined by the larger value and that defined by the smaller value.<br>
!>              That is, `delta` is the amount that the two input values need to be different in order to get a different result from the two weak and strong scaling methods.<br>
!>              Without loss of generality, suppose `x` is the larger value and that both `x` and `y` are positive. Therefore,<br>
!>              \code{.F90}
!>                  delta = reltol * (x - y)
!>                  delta / reltol = (x - y)
!>              \endcode
!>              The largest absolute difference that would pass the test, `(x - y)`, equals the tolerance times the larger value.<br>
!>              \code{.F90}
!>                  (x - y) = reltol * x
!>              \endcode
!>              Substituting into the expression for delta, one gets,<br>
!>              \code{.F90}
!>                  delta / reltol = reltol * x
!>                  delta = reltol**2 * x
!>              \endcode
!>              For example, if `x = 10; y = 9; reltol = 0.1`, then,<br>
!>              <ul>
!>                  <li>    The maximum tolerance is `reltol * x == 0.1 * 10 == 1.0`.<br>
!>                  <li>    The minimum tolerance is `reltol * y == 0.1 * 9.0 == 0.9`.<br>
!>              </ul>
!>              Therefore, `delta = (1.0 - 0.9) = 0.1` or `reltol**2 * x = 0.1**2 * 10 = .1`.<br>
!>              The absolute difference between the maximum and minimum tolerance tests in this case could be therefore, substantial.<br>
!>              However, the primary use case for the generic interfaces of this module is testing the results of computations.<br>
!>              In such cases, a relative tolerance is likely to be selected of much smaller magnitudes.<br>
!>              For example, a relative tolerance of `1e-8` is about half the precision available in a `real64`.<br>
!>              Therefore, the difference between the two tests would be `1e-8**2 * x or 1e-16 * x`, which is close to the limit of precision of `real64`.<br>
!>              As such, a relative tolerance of `reltol = 1e-9` (or smaller) would eliminate the difference between the two tests down to the limits of a `real64` precision.<br>
!>              That is, *each of the four methods will yield exactly the same results for all `x` and `y` values*.<br>
!>              <b>The symmetry of Relative Proximity</b><br>
!>              A relative comparison can be either **symmetric** or **asymmetric**.<br>
!>              <ul>
!>                  <li>    A *symmetric* relative comparison `isClose(x, y, reltol, abstol)` always yields the same result as `isClose(y, x)`.<br>
!>                  <li>    An *asymmetric* relative comparison is one uses only one of the values (as with the *Reference* comparison above).<br>
!>                          For an asymmetric comparison , `isClose(x, y, reltol, abstol)` is not necessarily the same as `isClose(y,x)`.<br>
!>              </ul>
!>              The suitability of symmetric vs. asymmetric methods depends on the question is being asked.<br>
!>              Symmetric tests are the most appropriate when,<br>
!>              <ul>
!>                  <li>    the question is: *Are the two numbers close to each other?* or,<br>
!>                  <li>    there is no obvious ordering between the two numbers.<br>
!>              </ul>
!>              Asymmetric tests are the most appropriate when,<br>
!>              <ul>
!>                  <li>    the question is: *Is the computed value `y` within `x` percentage of this known value `x`?* or,<br>
!>                  <li>    there is a **reference value** against which another value must be compared.<br>
!>              </ul>
!>              The symmetric approach provides an appealing consistency.<br>
!>              It mirrors the symmetry of equality, and is less confusing.<br>
!>              A symmetric test also eliminates the need for thinking about the order of the two input arguments `(x, b)`.<br>
!>              Conversely, when one needs to know that a value `y` is within a particular range of a known value `x`, then the following test can be implemented directly,<br>
!>              \code{.F90}
!>                  if (abs(x - y) <= reltol * x) then
!>                      ...
!>                  end if
!>              \endcode
!>              <b>Which symmetric comparison test is the most appropriate?</b><br>
!>              There are three possible symmetric tests as mentioned above:
!>              <ul>
!>                  <li>    <b>Weak Scaling</b>     : The maximum absolute value of the two.<br>
!>                          This approach yields the same result as the *Strong scaling* in most practical cases for small relative tolerance values.<br>
!>                  <li>    <b>Strong Scaling</b>   : The minimum absolute value of the two.<br>
!>                          This approach yields the same result as the *Strong scaling* in most practical cases for small relative tolerance values.<br>
!>                  <li>    <b>mean Scaling</b>     : The absolute value of the arithmetic mean of the two.<br>
!>                          This approach can potentially lead to the following **runtime complications**:<br>
!>                          <ul>
!>                              <li>    An arithmetic **overflow** can occur if the two values are large and are added together before dividing by two.<br>
!>                              <li>    An arithmetic **underflow** can occur if the two values are tiny and are divided by two before being added together.<br>
!>                          </ul>
!>              </ul>
!>              Despite the similarities of the *Weak* and *Strong* Scaling methods, the *Weak Scaling* approach provides a more useful result for very large relative tolerances.<br>
!>              This happens when one needs to test if two fairly disparate values are within a particular range of each other.<br>
!>              For example: *Is `x` within `200%` (`reltol = 2.0`) of `y`?<br>
!>              The weak testing would use the larger (non-zero) value for the test, and thus return `.true.` if one value is zero.<br>
!>              However, the Strong Scaling test would never indicate that two values are within that range of each other if one of them is zero.<br>
!>              For better illustration consider the following question: *Is `0` within `200%` of `10`?*<br>
!>              Note that `200%` of `10` is `20`. Therefore, the range within `200%` of `10` is `-10` to `+30`.<br>
!>              Therefore, zero falls within the range and the *Weak Scaling* test will return `.true.`.<br>
!>  </ol>
!>  <b>The logic behind the specific choices of absolute and relative tolerances.</b><br>
!>  <ul>
!>      <li>    The **relative tolerance** required for two values to be considered *close* is entirely use-case dependent.<br>
!>              Nevertheless, the relative tolerance needs to be greater than `epsilon(x)`.<br>
!>              This default tolerance is rather conservative knowing that the `reltol = sqrt(epsilon(x))` is typically the
!>              largest relative tolerance for which the various possible testing methods above will yield the similar results.<br>
!>              Also, good numerical algorithms are generally expected to not lose more than about half of available digits of accuracy due to errors.<br>
!>      <li>    The **absolute tolerance** value is primarily used for comparing to zero.<br>
!>              The absolute tolerance required to determine if a value is *close* to zero is entirely use-case dependent.<br>
!>              There is also essentially no bounds to the useful range. Thus a default of `tiny(x)` seems most appropriate.<br>
!>              This way, if a comparison with zero is to be made, the test is guaranteed to fail the first time, prompting the user to select an appropriate value subsequently.<br>
!>  </ul>
!>
!>  \see
!>  [math.isclose()](https://docs.python.org/3/library/math.html#math.isclose)<br>
!>  [numpy.isclose()](https://numpy.org/doc/stable/reference/generated/numpy.isclose.html)<br>
!>  [Floating-point comparison algorithms](https://www.boost.org/doc/libs/1_34_0/libs/test/doc/components/test_tools/floating_point_comparison.html)<br>
!>  [Python PEP 485](https://peps.python.org/pep-0485/#proposed-implementation)<br>
!>
!>  \test
!>  [test_pm_mathCompare](@ref test_pm_mathCompare)
!>
!>  \finmain
!>
!>  \author
!>  \FatemehBagheri, Thursday 12:45 AM, August 20, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_mathCompare

    use pm_kind, only: SK, IK, LK
    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathCompare"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is an empty derived type that is exclusively used to differentiate the procedures within the generic interface [isClose()](@ref pm_mathCompare::isClose).<br>
    !>
    !>  \details
    !>  Passing an object of this derived type to the procedures of the generic interface [isClose()](@ref pm_mathCompare::isClose)
    !>  is equivalent to requesting the **Reference Scaling** method of testing the proximity of two numbers.<br>
    !>
    !>  See the documentation of [pm_mathCompare](@ref pm_mathCompare) for extensive details of this and other possible testing modes.<br>
    !>  See the documentation of [isClose()](@ref pm_mathCompare::isClose) for example usage.<br>
    !>
    !>  \interface{reference_type}
    !>  \code{.F90}
    !>
    !>      use pm_mathCompare, only: reference_type
    !>      type(reference_type), parameter :: REFERENCE
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [isClose()](@ref pm_mathCompare::isClose)<br>
    !>
    !>  \finmain{reference_type}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 08:49 PM, August 10, 2021, Dallas, TX
    type :: reference_type
    end type

    !>  \brief
    !>  This is a convenience constant object of type [reference_type()](@ref pm_mathCompare::reference_type) that is exclusively provided
    !>  to distinguish and facilitate the use of the procedures under the generic interface [isClose()](@ref pm_mathCompare::isClose).<br>
    !>
    !>  \details
    !>  Passing this constant object to the procedures of the generic interface [isClose()](@ref pm_mathCompare::isClose)
    !>  is equivalent to requesting the **Reference Scaling** method of testing the proximity of two numbers.<br>
    !>
    !>  See the documentation of [pm_mathCompare](@ref pm_mathCompare) for extensive details of this and other possible testing modes.<br>
    !>  See the documentation of [isClose()](@ref pm_mathCompare::isClose) for example usage.<br>
    !>
    !>  \interface{REFERENCE}
    !>  \code{.F90}
    !>
    !>      use pm_mathCompare, only: REFERENCE
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [isClose()](@ref pm_mathCompare::isClose)<br>
    !>
    !>  \finmain{REFERENCE}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 08:49 PM, August 10, 2021, Dallas, TX
    type(reference_type), parameter :: REFERENCE = reference_type()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is an empty derived type that is exclusively used to differentiate the procedures within the generic interface [isClose()](@ref pm_mathCompare::isClose).<br>
    !>
    !>  \details
    !>  Passing an object of this derived type to the procedures of the generic interface [isClose()](@ref pm_mathCompare::isClose)
    !>  is equivalent to requesting the **Strong Scaling** method of testing the proximity of two numbers.<br>
    !>
    !>  See the documentation of [pm_mathCompare](@ref pm_mathCompare) for extensive details of this and other possible testing modes.<br>
    !>  See the documentation of [isClose()](@ref pm_mathCompare::isClose) for example usage.<br>
    !>
    !>  \interface{strong_type}
    !>  \code{.F90}
    !>
    !>      use pm_mathCompare, only: strong_type
    !>      type(strong_type), parameter :: STRONG
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [isClose()](@ref pm_mathCompare::isClose)<br>
    !>
    !>  \finmain{strong_type}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 08:49 PM, August 10, 2021, Dallas, TX
    type :: strong_type
    end type

    !>  \brief
    !>  This is a convenience constant object of type [strong_type()](@ref pm_mathCompare::strong_type) that is exclusively provided
    !>  to distinguish and facilitate the use of the procedures under the generic interface [isClose()](@ref pm_mathCompare::isClose).<br>
    !>
    !>  \details
    !>  Passing this constant object to the procedures of the generic interface [isClose()](@ref pm_mathCompare::isClose)
    !>  is equivalent to requesting the **Strong Scaling** method of testing the proximity of two numbers.<br>
    !>
    !>  See the documentation of [pm_mathCompare](@ref pm_mathCompare) for extensive details of this and other possible testing modes.<br>
    !>  See the documentation of [isClose()](@ref pm_mathCompare::isClose) for example usage.<br>
    !>
    !>  \interface{STRONG}
    !>  \code{.F90}
    !>
    !>      use pm_mathCompare, only: STRONG
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [isClose()](@ref pm_mathCompare::isClose)<br>
    !>
    !>  \finmain{STRONG}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 08:49 PM, August 10, 2021, Dallas, TX
    type(strong_type), parameter :: STRONG = strong_type()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is an empty derived type that is exclusively used to differentiate the procedures within the generic interface [isClose()](@ref pm_mathCompare::isClose).<br>
    !>
    !>  \details
    !>  Passing an object of this derived type to the procedures of the generic interface [isClose()](@ref pm_mathCompare::isClose)
    !>  is equivalent to requesting the **Weak Scaling** method of testing the proximity of two numbers.<br>
    !>
    !>  See the documentation of [pm_mathCompare](@ref pm_mathCompare) for extensive details of this and other possible testing modes.<br>
    !>  See the documentation of [isClose()](@ref pm_mathCompare::isClose) for example usage.<br>
    !>
    !>  \interface{weak_type}
    !>  \code{.F90}
    !>
    !>      use pm_mathCompare, only: weak_type
    !>      type(weak_type), parameter :: WEAK
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [isClose()](@ref pm_mathCompare::isClose)<br>
    !>
    !>  \finmain{weak_type}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 08:49 PM, August 10, 2021, Dallas, TX
    type :: weak_type
    end type

    !>  \brief
    !>  This is a convenience constant object of type [weak_type()](@ref pm_mathCompare::weak_type) that is exclusively provided
    !>  to distinguish and facilitate the use of the procedures under the generic interface [isClose()](@ref pm_mathCompare::isClose).<br>
    !>
    !>  \details
    !>  Passing this constant object to the procedures of the generic interface [isClose()](@ref pm_mathCompare::isClose)
    !>  is equivalent to requesting the **Weak Scaling** method of testing the proximity of two numbers.<br>
    !>
    !>  See the documentation of [pm_mathCompare](@ref pm_mathCompare) for extensive details of this and other possible testing modes.<br>
    !>  See the documentation of [isClose()](@ref pm_mathCompare::isClose) for example usage.<br>
    !>
    !>  \interface{WEAK}
    !>  \code{.F90}
    !>
    !>      use pm_mathCompare, only: WEAK
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [isClose()](@ref pm_mathCompare::isClose)<br>
    !>
    !>  \finmain{WEAK}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 08:49 PM, August 10, 2021, Dallas, TX
    type(weak_type), parameter :: WEAK = weak_type()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is an empty derived type that is exclusively used to differentiate the procedures within the generic interface [isClose()](@ref pm_mathCompare::isClose).<br>
    !>
    !>  \details
    !>  Passing an object of this derived type to the procedures of the generic interface [isClose()](@ref pm_mathCompare::isClose)
    !>  is equivalent to requesting the **mean Scaling** method of testing the proximity of two numbers.<br>
    !>
    !>  See the documentation of [pm_mathCompare](@ref pm_mathCompare) for extensive details of this and other possible testing modes.<br>
    !>  See the documentation of [isClose()](@ref pm_mathCompare::isClose) for example usage.<br>
    !>
    !>  \interface{mean_type}
    !>  \code{.F90}
    !>
    !>      use pm_mathCompare, only: mean_type
    !>      type(mean_type), parameter :: MEAN
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [isClose()](@ref pm_mathCompare::isClose)<br>
    !>
    !>  \finmain{mean_type}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 08:49 PM, August 10, 2021, Dallas, TX
    type :: mean_type
    end type

    !>  \brief
    !>  This is a convenience constant object of type [mean_type()](@ref pm_mathCompare::mean_type) that is exclusively provided
    !>  to distinguish and facilitate the use of the procedures under the generic interface [isClose()](@ref pm_mathCompare::isClose).<br>
    !>
    !>  \details
    !>  Passing this constant object to the procedures of the generic interface [isClose()](@ref pm_mathCompare::isClose)
    !>  is equivalent to requesting the **mean Scaling** method of testing the proximity of two numbers.<br>
    !>
    !>  See the documentation of [pm_mathCompare](@ref pm_mathCompare) for extensive details of this and other possible testing modes.<br>
    !>  See the documentation of [isClose()](@ref pm_mathCompare::isClose) for example usage.<br>
    !>
    !>  \interface{MEAN}
    !>  \code{.F90}
    !>
    !>      use pm_mathCompare, only: MEAN
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [isClose()](@ref pm_mathCompare::isClose)<br>
    !>
    !>  \finmain{MEAN}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 08:49 PM, August 10, 2021, Dallas, TX
    type(mean_type), parameter :: MEAN = mean_type()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the two input values are sufficiently close to each other within the specified tolerances.<br>
    !>
    !>  See the documentation of [pm_mathCompare](@ref pm_mathCompare) for extensive details of possible testing modes.<br>
    !>
    !>  \param[in]  x       :   The input scalar, or array of the same rank as other array-like arguments, of,
    !>                          <ul>
    !>                              <li>    type `complex` of kind \CKALL or, <br>
    !>                              <li>    type `real` of kind \RKALL,
    !>                          </ul>
    !>                          whose proximity to the input `y` will be tested based on the specified method and tolerances.<br>
    !>  \param[in]  y       :   The input scalar, or array of the same rank as other array-like arguments, of the same type and kind as the input `x`,
    !>                          whose proximity to the input `x` will be tested based on the specified method and tolerances.<br>
    !>  \param[in]  method  :   The input scalar, or array of the same rank as other array-like arguments, of,<br>
    !>                          <ul>
    !>                              <li>    type [reference_type](@ref pm_mathCompare::reference_type) or, <br>
    !>                              <li>    type [strong_type](@ref pm_mathCompare::strong_type) or, <br>
    !>                              <li>    type [weak_type](@ref pm_mathCompare::weak_type) or, <br>
    !>                              <li>    type [mean_type](@ref pm_mathCompare::mean_type), <br>
    !>                          </ul>
    !>                          whose presence is solely used to determine the scaling of the relative tolerance for the proximity test.<br>
    !>                          If `method = reference_type()`, then the input `x` will be taken as the reference value.<br>
    !>                          (**optional**, default = [WEAK](@ref pm_mathCompare::WEAK))
    !>  \param[in]  reltol  :   The input positive scalar, or array of the same rank as other array-like arguments, of type `real` of the same kind as the input argument `x`,
    !>                          representing the relative tolerance to be used in the proximity test.<br>
    !>                          (**optional**, default = `epsilon(x)`)
    !>  \param[in]  abstol  :   The input positive scalar, or array of the same rank as other array-like arguments, of type `real` of the same kind as the input argument `x`,
    !>                          representing the absolute tolerance to be used in the proximity test.<br>
    !>                          (**optional**, default = `tiny(x)`)
    !>
    !>  \return
    !>  `close`             :   The output scalar, or array of the same rank as other array-like arguments, of type `logical` of default kind \LK.<br>
    !>                          It is `.true.` <b>if and only if</b> the two input values `x` and `y` is close to each other according to the specified testing method and tolerances.<br>
    !>
    !>  \interface{isClose}
    !>  \code{.F90}
    !>
    !>      use pm_mathCompare, only: isClose
    !>      use pm_kind, only: LK
    !>      logical(LK) :: close
    !>
    !>      close = isClose(x, y, reltol = reltol, abstol = abstol)
    !>      close = isClose(x, y, method, reltol = reltol, abstol = abstol)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0. <= abstol` must hold for the corresponding arguments.<br>
    !>  The condition `0. <= reltol` must hold for the corresponding arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  A `-Inf` value is considered **not close** to a `+Inf` in all circumstances.<br>
    !>  Two `+Inf` values or two `-Inf` values are considered **close** to each other in all circumstances.<br>
    !>  A `NaN` value is considered **not close** to any other value including another `NaN` in all circumstances.<br>
    !>
    !>  \see
    !>  [pm_test](@ref pm_test)<br>
    !>
    !>  \example{isClose}
    !>  \include{lineno} example/pm_mathCompare/isClose/main.F90
    !>  \compilef{isClose}
    !>  \output{isClose}
    !>  \include{lineno} example/pm_mathCompare/isClose/main.out.F90
    !>
    !>  \test
    !>  [test_pm_mathCompare](@ref test_pm_mathCompare)
    !>
    !>  \todo
    !>  The implementations of the procedures of this generic interface can be improved
    !>  to ensure robustness against possible rare cases of overflows and underflow as
    !>  discussed in the documentation of [pm_mathCompare](@ref pm_mathCompare).<br>
    !>
    !>  \finmain{isClose}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 08:49 PM, August 10, 2021, Dallas, TX
    interface isClose

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE elemental module function isCloseDefault_CK5(x, y, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseDefault_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in)            :: x, y
        real(CKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if CK4_ENABLED
    PURE elemental module function isCloseDefault_CK4(x, y, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseDefault_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in)            :: x, y
        real(CKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if CK3_ENABLED
    PURE elemental module function isCloseDefault_CK3(x, y, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseDefault_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in)            :: x, y
        real(CKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if CK2_ENABLED
    PURE elemental module function isCloseDefault_CK2(x, y, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseDefault_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in)            :: x, y
        real(CKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if CK1_ENABLED
    PURE elemental module function isCloseDefault_CK1(x, y, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseDefault_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in)            :: x, y
        real(CKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function isCloseDefault_RK5(x, y, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseDefault_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in)            :: x, y
        real(RKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function isCloseDefault_RK4(x, y, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseDefault_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in)            :: x, y
        real(RKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function isCloseDefault_RK3(x, y, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseDefault_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in)            :: x, y
        real(RKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function isCloseDefault_RK2(x, y, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseDefault_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in)            :: x, y
        real(RKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function isCloseDefault_RK1(x, y, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseDefault_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in)            :: x, y
        real(RKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE elemental module function isCloseReference_CK5(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseReference_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in)            :: x, y
        type(reference_type), intent(in)            :: method
        real(CKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if CK4_ENABLED
    PURE elemental module function isCloseReference_CK4(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseReference_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in)            :: x, y
        type(reference_type), intent(in)            :: method
        real(CKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if CK3_ENABLED
    PURE elemental module function isCloseReference_CK3(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseReference_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in)            :: x, y
        type(reference_type), intent(in)            :: method
        real(CKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if CK2_ENABLED
    PURE elemental module function isCloseReference_CK2(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseReference_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in)            :: x, y
        type(reference_type), intent(in)            :: method
        real(CKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if CK1_ENABLED
    PURE elemental module function isCloseReference_CK1(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseReference_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in)            :: x, y
        type(reference_type), intent(in)            :: method
        real(CKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function isCloseReference_RK5(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseReference_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in)            :: x, y
        type(reference_type), intent(in)            :: method
        real(RKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function isCloseReference_RK4(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseReference_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in)            :: x, y
        type(reference_type), intent(in)            :: method
        real(RKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function isCloseReference_RK3(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseReference_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in)            :: x, y
        type(reference_type), intent(in)            :: method
        real(RKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function isCloseReference_RK2(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseReference_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in)            :: x, y
        type(reference_type), intent(in)            :: method
        real(RKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function isCloseReference_RK1(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseReference_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in)            :: x, y
        type(reference_type), intent(in)            :: method
        real(RKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE elemental module function isCloseStrong_CK5(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseStrong_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in)            :: x, y
        type(strong_type)   , intent(in)            :: method
        real(CKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if CK4_ENABLED
    PURE elemental module function isCloseStrong_CK4(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseStrong_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in)            :: x, y
        type(strong_type)   , intent(in)            :: method
        real(CKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if CK3_ENABLED
    PURE elemental module function isCloseStrong_CK3(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseStrong_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in)            :: x, y
        type(strong_type)   , intent(in)            :: method
        real(CKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if CK2_ENABLED
    PURE elemental module function isCloseStrong_CK2(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseStrong_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in)            :: x, y
        type(strong_type)   , intent(in)            :: method
        real(CKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if CK1_ENABLED
    PURE elemental module function isCloseStrong_CK1(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseStrong_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in)            :: x, y
        type(strong_type)   , intent(in)            :: method
        real(CKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function isCloseStrong_RK5(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseStrong_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in)            :: x, y
        type(strong_type)   , intent(in)            :: method
        real(RKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function isCloseStrong_RK4(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseStrong_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in)            :: x, y
        type(strong_type)   , intent(in)            :: method
        real(RKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function isCloseStrong_RK3(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseStrong_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in)            :: x, y
        type(strong_type)   , intent(in)            :: method
        real(RKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function isCloseStrong_RK2(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseStrong_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in)            :: x, y
        type(strong_type)   , intent(in)            :: method
        real(RKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function isCloseStrong_RK1(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseStrong_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in)            :: x, y
        type(strong_type)   , intent(in)            :: method
        real(RKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE elemental module function isCloseWeak_CK5(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseWeak_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in)            :: x, y
        type(weak_type)     , intent(in)            :: method
        real(CKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if CK4_ENABLED
    PURE elemental module function isCloseWeak_CK4(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseWeak_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in)            :: x, y
        type(weak_type)     , intent(in)            :: method
        real(CKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if CK3_ENABLED
    PURE elemental module function isCloseWeak_CK3(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseWeak_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in)            :: x, y
        type(weak_type)     , intent(in)            :: method
        real(CKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if CK2_ENABLED
    PURE elemental module function isCloseWeak_CK2(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseWeak_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in)            :: x, y
        type(weak_type)     , intent(in)            :: method
        real(CKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if CK1_ENABLED
    PURE elemental module function isCloseWeak_CK1(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseWeak_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in)            :: x, y
        type(weak_type)     , intent(in)            :: method
        real(CKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function isCloseWeak_RK5(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseWeak_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in)            :: x, y
        type(weak_type)     , intent(in)            :: method
        real(RKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function isCloseWeak_RK4(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseWeak_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in)            :: x, y
        type(weak_type)     , intent(in)            :: method
        real(RKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function isCloseWeak_RK3(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseWeak_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in)            :: x, y
        type(weak_type)     , intent(in)            :: method
        real(RKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function isCloseWeak_RK2(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseWeak_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in)            :: x, y
        type(weak_type)     , intent(in)            :: method
        real(RKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function isCloseWeak_RK1(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseWeak_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in)            :: x, y
        type(weak_type)     , intent(in)            :: method
        real(RKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE elemental module function isCloseMean_CK5(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseMean_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in)            :: x, y
        type(mean_type)     , intent(in)            :: method
        real(CKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if CK4_ENABLED
    PURE elemental module function isCloseMean_CK4(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseMean_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)        , intent(in)            :: x, y
        type(mean_type)     , intent(in)            :: method
        real(CKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if CK3_ENABLED
    PURE elemental module function isCloseMean_CK3(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseMean_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)        , intent(in)            :: x, y
        type(mean_type)     , intent(in)            :: method
        real(CKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if CK2_ENABLED
    PURE elemental module function isCloseMean_CK2(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseMean_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)        , intent(in)            :: x, y
        type(mean_type)     , intent(in)            :: method
        real(CKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if CK1_ENABLED
    PURE elemental module function isCloseMean_CK1(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseMean_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)        , intent(in)            :: x, y
        type(mean_type)     , intent(in)            :: method
        real(CKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function isCloseMean_RK5(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseMean_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)           , intent(in)            :: x, y
        type(mean_type)     , intent(in)            :: method
        real(RKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function isCloseMean_RK4(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseMean_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)           , intent(in)            :: x, y
        type(mean_type)     , intent(in)            :: method
        real(RKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function isCloseMean_RK3(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseMean_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)           , intent(in)            :: x, y
        type(mean_type)     , intent(in)            :: method
        real(RKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function isCloseMean_RK2(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseMean_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)           , intent(in)            :: x, y
        type(mean_type)     , intent(in)            :: method
        real(RKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function isCloseMean_RK1(x, y, method, reltol, abstol) result(close)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCloseMean_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)           , intent(in)            :: x, y
        type(mean_type)     , intent(in)            :: method
        real(RKC)           , intent(in), optional  :: reltol, abstol
        logical(LK)                                 :: close
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface isClose

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_mathCompare ! LCOV_EXCL_LINE