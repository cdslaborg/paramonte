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
!>  This module contains relevant mathematical constants.
!>
!>  \note
!>  The constants of this module are saved with the highest available `real` precision kind.<br>
!>  To use the constants at expressions involving lower-precision `real` kinds, simply convert the numbers to the desired kind
!>  via the Fortran intrinsic `real(x, kind = RKG)` where `RKG` refers to the target kind parameter used in the expression.
!>
!>  \final
!>
!>  \todo
!>  \plow
!>  Th contents of this module can be expanded to include more mathematical constants.<br>
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_mathConst

    use pm_kind, only: SK, RKB

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_mathConst"

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! \pi-related constants.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    real(RKB)   , parameter :: PI = acos(-1._RKB)                                       !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the irrational number \f$\pi\f$.
    real(RKB)   , parameter :: TWO_PI = 2 * PI                                          !<  \public The scalar `real` constant of kind with highest available precision \RKB representing twice the irrational number \f$\pi\f$.
    real(RKB)   , parameter :: HALF_PI = .5_RKB * PI                                    !<  \public The scalar `real` constant of kind with highest available precision \RKB representing half the irrational number \f$\pi\f$.
    real(RKB)   , parameter :: INVERSE_PI = 1._RKB / PI                                 !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the inverse of the irrational number \f$\pi\f$.
    real(RKB)   , parameter :: QUARTER_PI = .25_RKB * PI                                !<  \public The scalar `real` constant of kind with highest available precision \RKB representing a quarter of the irrational number \f$\pi\f$.
    real(RKB)   , parameter :: LOG_PI = log(PI)                                         !<  \public The scalar `real` constant of kind with highest available precision \RKB representing \f$\log(\pi)\f$.
    real(RKB)   , parameter :: SQRT_PI = sqrt(PI)                                       !<  \public The scalar `real` constant of kind with highest available precision \RKB representing \f$\sqrt{\pi}\f$.
    real(RKB)   , parameter :: LOG_TWO_PI = log(TWO_PI)                                 !<  \public The scalar `real` constant of kind with highest available precision \RKB representing \f$\log(2\pi)\f$.
    real(RKB)   , parameter :: SQRT_TWO_PI = sqrt(TWO_PI)                               !<  \public The scalar `real` constant of kind with highest available precision \RKB representing \f$\sqrt{2\pi}\f$.
    real(RKB)   , parameter :: SQRT_HALF_PI = sqrt(HALF_PI)                             !<  \public The scalar `real` constant of kind with highest available precision \RKB representing \f$\sqrt{\frac{\pi}{2}}\f$.
    real(RKB)   , parameter :: INVERSE_SQRT_PI = sqrt(INVERSE_PI)                       !<  \public The scalar `real` constant of kind with highest available precision \RKB representing \f$\frac{1}{\sqrt{\pi}}\f$.
    real(RKB)   , parameter :: INVERSE_SQRT_TWO_PI = 1._RKB / SQRT_TWO_PI               !<  \public The scalar `real` constant of kind with highest available precision \RKB representing \f$\frac{1}{\sqrt{2\pi}}\f$.
    real(RKB)   , parameter :: LOG_INVERSE_SQRT_TWO_PI = log(INVERSE_SQRT_TWO_PI)       !<  \public The scalar `real` constant of kind with highest available precision \RKB representing \f$\log\left(\frac{1}{\sqrt{2\pi}}\right)\f$, frequently appearing in distributions (e.g., [Normal](@ref pm_distNorm)).

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Common numeric constants.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    real(RKB)    , parameter :: NAPIER = exp(1._RKB)                                    !<  \public The scalar `real` constant of kind with highest available precision \RKB representing the Napier constant (a.k.a. Euler number) \f$e = \exp(1)\f$.
    real(RKB)    , parameter :: LOG_TWO = log(2._RKB)                                   !<  \public The scalar `real` constant of kind with highest available precision \RKB representing \f$\log(2)\f$.
    real(RKB)    , parameter :: LOG_TEN = log(1.e1_RKB)                                 !<  \public The scalar `real` constant of kind with highest available precision \RKB representing \f$\log(10)\f$.
    real(RKB)    , parameter :: LOG_HALF = log(0.5_RKB)                                 !<  \public The scalar `real` constant of kind with highest available precision \RKB representing \f$\log\left(\frac{1}{2}\right)\f$.
    real(RKB)    , parameter :: SQRT_TWO = sqrt(2._RKB)                                 !<  \public The scalar `real` constant of kind with highest available precision \RKB representing \f$\sqrt{2}\f$.
    real(RKB)    , parameter :: LOG10_NAPIER = log10(NAPIER)                            !<  \public The scalar `real` constant of kind with highest available precision \RKB representing \f$\log_{10}(e)\f$.
    real(RKB)    , parameter :: INVERSE_LOG_TWO = 1._RKB / LOG_TWO                      !<  \public The scalar `real` constant of kind with highest available precision \RKB representing \f$\frac{1}{\log(2)}\f$.
    real(RKB)    , parameter :: INVERSE_SQRT_TWO = 1._RKB / SQRT_TWO                    !<  \public The scalar `real` constant of kind with highest available precision \RKB representing \f$\frac{1}{\sqrt{2}}\f$.

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the derived type [origin_type](@ref pm_mathConst::origin_type) representing the geometric origin of the coordinates.
    !>
    !>  \details
    !>  For example usage, see the corresponding PaaMonte generic interfaces that use this object.
    !>
    !>  The origin of a Euclidean space is a special point, usually denoted by the letter **O**, used as a fixed point of reference for the geometry of the surrounding space.<br>
    !>  In physical problems, the choice of origin is often arbitrary, meaning any choice of origin will ultimately give the same answer.<br>
    !>  This allows one to pick an origin point that makes the mathematics as simple as possible, often by taking advantage of some kind of geometric symmetry.<br>
    !>
    !>  In a Cartesian coordinate system, the origin is the point where the axes of the system intersect.<br>
    !>  The origin divides each of these axes into two halves, a positive and a negative semiaxis.<br>
    !>  Points can then be located with reference to the origin by giving their numerical coordinates—that is,
    !>  the positions of their projections along each axis, either in the positive or negative direction.<br>
    !>  **The coordinates of the origin are always all zero**, for example \f$(0,0)\f$ in two dimensions and \f$(0,0,0)\f$ in three.<br>
    !>
    !>  The origin of the complex plane can be referred as the point where real axis and imaginary axis intersect each other.<br>
    !>  In other words, **the origin in the complex plane is the complex number zero**.<br>
    !>
    !>  \see
    !>  [ORIGIN](@ref pm_mathConst::ORIGIN)<br>
    !>
    !>  \final{origin_type}
    !>
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    type :: origin_type; end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The scalar constant object of type [origin_type](@ref pm_mathConst::origin_type) representing the geometric origin of the coordinates.
    !>
    !>  \details
    !>  For example usage, see the corresponding PaaMonte generic interfaces that use this object.
    !>
    !>  The origin of a Euclidean space is a special point, usually denoted by the letter **O**, used as a fixed point of reference for the geometry of the surrounding space.<br>
    !>  In physical problems, the choice of origin is often arbitrary, meaning any choice of origin will ultimately give the same answer.<br>
    !>  This allows one to pick an origin point that makes the mathematics as simple as possible, often by taking advantage of some kind of geometric symmetry.<br>
    !>
    !>  In a Cartesian coordinate system, the origin is the point where the axes of the system intersect.<br>
    !>  The origin divides each of these axes into two halves, a positive and a negative semiaxis.<br>
    !>  Points can then be located with reference to the origin by giving their numerical coordinates—that is,
    !>  the positions of their projections along each axis, either in the positive or negative direction.<br>
    !>  **The coordinates of the origin are always all zero**, for example \f$(0,0)\f$ in two dimensions and \f$(0,0,0)\f$ in three.<br>
    !>
    !>  The origin of the complex plane can be referred as the point where real axis and imaginary axis intersect each other.<br>
    !>  In other words, **the origin in the complex plane is the complex number zero**.<br>
    !>
    !>  \see
    !>  [origin_type](@ref pm_mathConst::origin_type)<br>
    !>
    !>  \final{ORIGIN}
    !>
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    type(origin_type), parameter :: ORIGIN = origin_type()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The scalar `real` constant of kind with highest available precision \RKB representing the **Euler-Mascheroni** constant.
    !>
    !>  \details
    !>  The Euler constant (sometimes also called the Euler–Mascheroni constant) is a mathematical constant usually denoted by the lowercase Greek letter \f$\gamma\f$.
    !>  It is defined as the limiting difference between the harmonic series and the natural logarithm,
    !>  \f{equation}{
    !>  \begin{aligned}
    !>      \gamma
    !>      &= \lim _{n \to \infty} \left( -\log(n) + \sum _{k=1}^{n} \frac{1}{k} \right) \\
    !>      &= \int _{1}^{\infty} \left( -\frac{1}{x} + \frac{1}{\lfloor x \rfloor} \right) ~ dx ~.
    !>  \end{aligned}
    !>  \f}
    !>  Here, \f$\lfloor x \rfloor\f$ represents the Fortran intrinsic `floor()` function.
    !>
    !>  \image html pm_mathGamma@EULER_CONST.png "The area of the blue region converges to the Euler constant" width=300
    !>
    !>  <b>Applications</b><br>
    !>  The Euler constant appears in many areas of science, including,
    !>  <ol>
    !>      <li>    Expressions involving the exponential integral.
    !>      <li>    The Laplace transform of the natural logarithm.
    !>      <li>    The first term of the Laurent series expansion for the Riemann zeta function, where it is the first of the Stieltjes constants.
    !>      <li>    Calculations of the Digamma function.
    !>      <li>    A product formula for the Gamma function.
    !>      <li>    The asymptotic expansion of the Gamma function for small arguments.
    !>      <li>    An inequality for the Euler totient function.
    !>      <li>    The growth rate of the divisor function.
    !>      <li>    In dimensional regularization of Feynman diagrams in quantum field theory.
    !>      <li>    The calculation of the Meissel–Mertens constant.
    !>      <li>    The third of the Mertens theorems.
    !>      <li>    Solution of the second kind to the Bessel equation.
    !>      <li>    In the regularization/renormalization of the harmonic series as a finite value.
    !>      <li>    The mean of the Gumbel distribution.
    !>      <li>    The information entropy of the Weibull and Lévy distributions, and, implicitly, of the Chi-squared distribution for one or two degrees of freedom.
    !>      <li>    The answer to the Coupon Collector problem.
    !>      <li>    In some formulations of the Zipf law.
    !>      <li>    A definition of the cosine integral.
    !>      <li>    Lower bounds to a prime gap.
    !>      <li>    An upper bound on Shannon entropy in quantum information theory.
    !>      <li>    Fisher–Orr model for genetics of adaptation in evolutionary biology.
    !>  </ol>
    !>
    !>  \see
    !>  [The Euler Number](@ref pm_mathConst::NAPIER)<br>
    !>  [The Euler Constant](https://en.wikipedia.org/wiki/Euler%27s_constant)<br>
    !>
    !>  \todo
    !>  \pvhigh
    !>  \fixPrecLoss{this constant, the first `99`}
    !>
    !>  \final{EULER_CONST}
    !>
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    real(RKB)    , parameter :: EULER_CONST = 0.577215664901532860606512090082402431042159335939923598805767234884867726777664670936947063291746749_RKB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The scalar `real` constant of kind with highest available precision \RKB representing the **Apery** constant.
    !>
    !>  \details
    !>  In mathematics, Apéry's constant is the sum of the reciprocals of the positive cubes. That is, it is defined as the number
    !>  \f{equation}{
    !>      \begin{aligned}
    !>          \zeta(3)
    !>          &= \sum _{n=1}^{\infty} \frac{1}{n^{3}} \\
    !>          &= \lim _{n\to\infty} \left( {\frac{1}{1^{3}}} + {\frac {1}{2^{3}}} + \cdots + {\frac{1}{n^{3}}} \right) ~,
    !>      \end{aligned}
    !>  \f}
    !>  where \f$\zeta\f$ is the **Riemann zeta function**.<br>
    !>  It has an approximate value of,
    !>  \f{equation}{
    !>      \zeta(3) = 1.202056903159594285399738161511449990764986292\ldots ~.
    !>  \f}
    !>  The constant is named after **Roger Apéry**.<br>
    !>  It arises naturally in a number of physical problems, including in the second- and third-order terms of the electron's gyromagnetic ratio using **quantum electrodynamics**.<br>
    !>  It also arises in the analysis of random minimum spanning trees and in conjunction with the Gamma function when solving certain integrals involving exponential functions in a quotient.<br>
    !>  These appear occasionally in physics, for instance, when evaluating the two-dimensional case of the Debye model and the Stefan–Boltzmann law.<br>
    !>
    !>  \todo
    !>  \pvhigh
    !>  \fixPrecLoss{this constant, the first `109`}<br>
    !>
    !>  \final{APERY_CONST}
    !>
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    real(RKB)    , parameter :: APERY_CONST = 1.20205690315959428539973816151144999076498629234049888179227155534183820578631309018645587360933525814619915_RKB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The scalar `real` constant of kind with highest available precision \RKB representing the irrational **Prime** constant.
    !>
    !>  \details
    !>  The prime constant is the real number \f$\rho\f$ whose \f$n\f$th binary digit is \f$1\f$ if \f$n\f$ is prime and `0` if \f$n\f$ is composite or \f$1\f$.<br>
    !>  In other words, \f$\rho\f$ is the number whose binary expansion corresponds to the indicator function of the set of prime numbers.<br>
    !>  That is,
    !>  \f{equation}{
    !>      \rho = \sum_{p} \frac{1}{2^{p}} = \sum_{n = 1}^{\infty} \frac{\chi_{\mathbb{P}}(n)}{2^{n}} \rho = \sum_{{p}} \frac{1}{2^{p}} = \sum_{{n = 1}}^{\infty} \frac{\chi_{{{\mathbb{P}}}}(n)}{2^{n}} ~,
    !>  \f}
    !>  where \f$p\f$ indicates a prime and \f$\chi_{{{\mathbb{P}}}}\f$ is the characteristic function of the set \f$\mathbb{P}\f$ of prime numbers.<br>
    !>  The beginning of the decimal expansion of \f$\rho\f$ is: \f$\rho = 0.414682509851111660248109622\ldots\f$.<br>
    !>
    !>  \todo
    !>  \pvhigh
    !>  \fixPrecLoss{this constant, the first `105`}<br>
    !>
    !>  \final{PRIME_CONST}
    !>
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    real(RKB)    , parameter :: PRIME_CONST = .414682509851111660248109622154307708365774238137916977868245414488640960619357334196290048428475777939616_RKB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The scalar `real` constant of kind with highest available precision \RKB representing the **Golden Ratio** constant.
    !>
    !>  \details
    !>  In mathematics, two quantities are in the Golden Ratio if their ratio is the same as the ratio of their sum to the larger of the two quantities.<br>
    !>  Expressed algebraically, for quantities \f$a\f$ and \f$b\f$ with \f$0 < a < b\f$,
    !>  \f{equation}{
    !>      \frac{a + b}{a} = \frac{a}{b} = \phi ~,
    !>  \f}
    !>  where the Greek letter \f$\phi\f$ denotes the Golden Ratio.<br>
    !>  The constant \f$\phi\f$ satisfies the quadratic equation \f$\phi^{2} = \phi + 1\f$ and is an irrational number with a value of,
    !>  \f{equation}{
    !>      \phi = \frac{1 + \sqrt{5}}{2} = 1.618033988749\ldots ~.
    !>  \f}
    !>  The Golden Ratio was called the **extreme and mean ratio** by Euclid, and **the divine proportion** by Luca Pacioli, among other names.<br>
    !>  Mathematicians have studied the properties of the Golden Ratio since antiquity.<br>
    !>  It is the ratio of a diagonal of a regular pentagon to its side and thus appears in the construction of the dodecahedron and icosahedron.<br>
    !>  A **golden rectangle** — that is, a rectangle with an aspect ratio of \f$\phi\f$ — may be cut into a square and a smaller rectangle with the same aspect ratio.<br>
    !>  The Golden Ratio has been used to analyze the proportions of natural objects and artificial systems such as financial markets, in some cases based on dubious fits to data.<br>
    !>  <b>The Golden Ratio appears in some patterns in nature</b>, including the spiral arrangement of leaves and other parts of vegetation.<br>
    !>  Some 20th-century artists and architects, including Le Corbusier and Salvador Dalí, have proportioned their works to approximate the Golden Ratio, believing it to be aesthetically pleasing.<br>
    !>  These usages frequently appear in the form of a golden rectangle, as illustrated below.
    !>
    !>  \image html pm_mathGamma@GOLDEN_RATIO.png "A golden rectangle with long side `a` and short side `b` (shaded red, right) and a square with sides of length `a` (shaded blue, left) combine to form a similar golden rectangle with long side `a + b` and short side `a`." width=300
    !>
    !>  \see
    !>  [SILVER_RATIO](@ref pm_mathConst::SILVER_RATIO)<br>
    !>  [GOLDEN_RATIO](@ref pm_mathConst::GOLDEN_RATIO)<br>
    !>  [SUPER_GOLDEN_RATIO](@ref pm_mathConst::SUPER_GOLDEN_RATIO)<br>
    !>
    !>  \final{GOLDEN_RATIO}
    !>
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    real(RKB)    , parameter :: GOLDEN_RATIO = .5_RKB * (1._RKB + sqrt(5._RKB))

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The scalar `real` constant of kind with highest available precision \RKB representing the **Silver Ratio** constant.
    !>
    !>  \details
    !>  In mathematics, two quantities are in the **Silver Ratio** (or **silver mean**) if the ratio of the smaller of
    !>  those two quantities to the larger quantity is the same as the ratio of the larger quantity to the sum of the smaller quantity and twice the larger quantity.<br>
    !>  This defines the Silver Ratio as an irrational mathematical constant, whose value of one plus the square root of `2` is approximately `2.4142135623`.<br>
    !>  Its name is an allusion to the [Golden Ratio](@ref pm_mathConst::GOLDEN_RATIO).<br>
    !>  Analogously to the way the Golden Ratio is the limiting ratio of consecutive Fibonacci numbers,
    !>  the Silver Ratio is the limiting ratio of consecutive **Pell numbers**.<br>
    !>  The Silver Ratio is denoted by \f$\delta_S\f$, which can be expressed algebraically as,
    !>  \f{equation}{
    !>      \frac {2a+b}{a} = \frac {a}{b} \equiv \delta _{S} ~,
    !>  \f}
    !>  or equivalently,
    !>  \f{equation}{
    !>      2 +  \frac {b}{a} = {\frac {a}{b}} \equiv \delta _{S} ~.
    !>  \f}
    !>  The Silver Ratio can also be defined by the simple continued fraction \f$[2; 2, 2, 2, ...]\f$,
    !>  \f{equation}{
    !>      2 + {\cfrac{1}{2 + {\cfrac{1}{2 + {\cfrac{1}{2 + \ddots}}}}}} = \delta _{S} ~.
    !>  \f}
    !>  The convergents of this continued fraction are ratios of consecutive **Pell numbers**.<br>
    !>  These fractions provide accurate rational approximations of the Silver Ratio,
    !>  analogous to the approximation of the Golden Ratio by ratios of consecutive Fibonacci numbers.<br>
    !>  The **Silver Rectangle** is connected to the **regular octagon**.<br>
    !>  If a regular octagon is partitioned into two isosceles trapezoids and a rectangle, then the rectangle is a **Silver Rectangle** with an aspect ratio of \f$1:\delta_S\f$.<br>
    !>  The four sides of the trapezoids are in a ratio of \f$1:1:1:\delta_S\f$.<br>
    !>  If the edge length of a regular octagon is \f$t\f$, then the span of the octagon (the distance between opposite sides) is \f$\delta_{S}t\f$, and the area of the octagon is \f$2\delta_{S}t^2\f$.<br>
    !>
    !>  \image html pm_mathGamma@SILVER_RATIO.png "Silver Rectangle" width=300
    !>  <br>
    !>  \image html pm_mathGamma@SILVER_RATIO@octagon.png "Silver Ratio within the octagon." width=300
    !>
    !>  \see
    !>  [SILVER_RATIO](@ref pm_mathConst::SILVER_RATIO)<br>
    !>  [GOLDEN_RATIO](@ref pm_mathConst::GOLDEN_RATIO)<br>
    !>  [SUPER_GOLDEN_RATIO](@ref pm_mathConst::SUPER_GOLDEN_RATIO)<br>
    !>
    !>  \final{SILVER_RATIO}
    !>
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    real(RKB)    , parameter :: SILVER_RATIO = 1._RKB + sqrt(2._RKB)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The scalar `real` constant of kind with highest available precision \RKB representing the **Supergolden Ratio** constant.<br>
    !>
    !>  \details
    !>  In mathematics, two quantities are in the superGolden Ratio if the quotient of the larger number divided by the smaller one is equal to
    !>  \f{equation}{
    !>      \psi = {\frac{1 + {\sqrt[{3}]{\frac{29 + 3{\sqrt{93}}}{2}}} + {\sqrt[{3}]{\frac{29 - 3{\sqrt{93}}}{2}}}}{3}} ~,
    !>  \f}
    !>  which is the only real solution to the equation \f$x^{3} = x^{2} + 1\f$.<br>
    !>  It can also be represented using the hyperbolic cosine as,
    !>  \f{equation}{
    !>      \psi = \frac{2}{3} \cosh{\left( \tfrac{\cosh^{-1} \left(\frac{29}{2} \right)}{3} \right)} + \frac {1}{3} ~.
    !>  \f}
    !>  The decimal expansion of this number begins \f$1.465571231876768026656731\ldots\f$,
    !>  and the ratio is commonly represented by the Greek letter \f$\psi\f$.<br>
    !>
    !>  Many of the properties of the Supergolden Ratio are related to those of the [Golden Ratio](@ref pm_mathConst::GOLDEN_RATIO).<br>
    !>  For example, the \f$n\f$th item of the Narayana sequence is the number of ways to tile a \f$1\times n\f$ rectangle with \f$1\times 1\f$ and \f$1\times3\f$ tiles.<br>
    !>  Similarly, the \f$n\f$th term of the Fibonacci sequence is the number of ways to tile a \f$1\times n\f$ rectangle with \f$1\times1\f$ and \f$1\times2\f$ tiles.<br>
    !>  The Supergolden Ratio satisfies \f$\psi -1 = \psi^{-2}\f$, while the Golden Ratio satisfies \f$\phi -1 = \phi^{-1}\f$.<br>
    !>  In the Fibonacci Rabbit Problem, each pair breeds each cycle starting after two cycles, while in the Narayana Cow Problem, each pair breeds each cycle starting after three cycles.<br>
    !>  There is a **Supergolden Rectangle** that has the property that if a square is removed from one side, the remaining rectangle can be divided into two Supergolden Rectangles of opposite orientations.<br>
    !>
    !>  A **Supergolden Rectangle** is a rectangle whose side lengths are in the Supergolden Ratio, i.e.
    !>  the length of the longer side divided by the length of the shorter side is equal to the Supergolden Ratio: \f${\frac{1 + {\sqrt[{3}]{\frac{29 + 3{\sqrt{93}}}{2}}} + {\sqrt[{3}]{\frac{29-3{\sqrt{93}}}{2}}}}{3}}\f$.<br>
    !>  When a square with the same side length as the shorter side of the rectangle is removed from one side of the rectangle, the sides resulting rectangle will be in a \f$\psi^2:1\f$ ratio.<br>
    !>  This rectangle can be divided into rectangles with side-length ratios of \f$\psi:1\f$ and \f$1:\psi\f$, two Supergolden Ratios of perpendicular orientations.<br>
    !>  Their areas will be in a \f$\psi^2:1\f$ ratio.<br>
    !>  In addition, if the line that separates the two Supergolden Rectangles from each other is extended across the rest of the original rectangle
    !>  such that it – along with the side of the square that was removed from the original rectangle – divides the original rectangle into quadrants,
    !>  then the larger Supergolden Rectangle has the same area as the opposite quadrant,
    !>  its diagonal length is the length of the short side of the original rectangle divided by \f$\sqrt{\psi}\f$, the fourth quadrant is also a Supergolden Rectangle,
    !>  and its diagonal length is \f$\sqrt {\psi }\f$ times the length of the short side of the original rectangle.<br>
    !>
    !>  \image html pm_mathGamma@SUPER_GOLDEN_RATIO@triangle.png "A triangle with sides lengths of the Supergolden Ratio, its inverse, and one has an angle of exactly 120 degrees opposite the ratio length." width=300
    !>  <br>
    !>  \image html pm_mathGamma@SUPER_GOLDEN_RATIO@rectangle.png "This diagram shows the lengths of decreasing powers within a Supergolden Rectangle, and the pattern of intersecting right angles that appears as a result." width=300
    !>
    !>  \see
    !>  [SILVER_RATIO](@ref pm_mathConst::SILVER_RATIO)<br>
    !>  [GOLDEN_RATIO](@ref pm_mathConst::GOLDEN_RATIO)<br>
    !>  [SUPER_GOLDEN_RATIO](@ref pm_mathConst::SUPER_GOLDEN_RATIO)<br>
    !>
    !>  \final{SUPER_GOLDEN_RATIO}
    !>
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    real(RKB)    , parameter :: SUPER_GOLDEN_RATIO = (2._RKB * cosh(acosh(29._RKB / 2._RKB) / 3._RKB) + 1._RKB) / 3._RKB

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the indicator type for generating instances of objects that indicate the use of the
    !>  negative infinity \f$-\infty\f$ as an input argument to the generic interfaces of the ParaMonte library.<br>
    !>
    !>  \details
    !>  This is an empty derived type that exists solely for generating unique objects that are distinguishable
    !>  as input arguments to procedures under the generic interfaces of the ParaMonte library.<br>
    !>
    !>  \interface{ninf_type}
    !>  \code{.F90}
    !>
    !>      use pm_mathConst, only: ninf_type
    !>      type(ninf_type), parameter :: ninf
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [pm_quadPack](@ref pm_quadPack)<br>
    !>
    !>  \final{ninf_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: ninf_type
    end type

    !>  \brief
    !>  The scalar constant object of type [ninf_type](@ref pm_mathConst::ninf_type) that indicates the use of the
    !>  negative infinity \f$-\infty\f$ as an input argument to the generic interfaces of the ParaMonte library.<br>
    !>
    !>  \details
    !>  See the documentations of the generic interfaces that use this constant for example usage.<br>
    !>
    !>  \see
    !>  [pm_quadPack](@ref pm_quadPack)<br>
    !>
    !>  \final{ninf}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(ninf_type) , parameter :: ninf = ninf_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: ninf
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is the indicator type for generating instances of objects that indicate the use of the
    !>  positive infinity \f$+\infty\f$ as an input argument to the generic interfaces of the ParaMonte library.<br>
    !>
    !>  \details
    !>  This is an empty derived type that exists solely for generating unique objects that are distinguishable
    !>  as input arguments to procedures under the generic interfaces of the ParaMonte library.<br>
    !>
    !>  \interface{pinf_type}
    !>  \code{.F90}
    !>
    !>      use pm_mathConst, only: pinf_type
    !>      type(pinf_type), parameter :: pinf
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [pm_quadPack](@ref pm_quadPack)<br>
    !>
    !>  \final{pinf_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: pinf_type
    end type

    !>  \brief
    !>  The scalar constant object of type [pinf_type](@ref pm_mathConst::pinf_type) that indicates the use of the
    !>  positive infinity \f$+\infty\f$ as an input argument to the generic interfaces of the ParaMonte library.<br>
    !>
    !>  \details
    !>  See the documentations of the generic interfaces that use this constant for example usage.<br>
    !>
    !>  \see
    !>  [pm_quadPack](@ref pm_quadPack)<br>
    !>
    !>  \final{pinf}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(pinf_type) , parameter :: pinf = pinf_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: pinf
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_mathConst ! LCOV_EXCL_LINE