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
!>  This module contains classes and procedures for computing the roots of
!>  one-dimensional continuous mathematical functions using various root-finding methods.<br>
!>
!>  \details
!>
!>  In mathematics and computing, a **root-finding** algorithm is an algorithm for finding **zeros**, also called **roots**, of continuous functions.<br>
!>  A zero of a function \f$f\f$, from the real numbers to real numbers or from the complex numbers to the complex numbers, is a number \f$x\f$ such that \f$f(x) = 0\f$.<br>
!>  As, generally, the zeros of a function cannot be computed exactly nor expressed in closed form, root-finding algorithms provide approximations to zeros
!>  expressed either as floating-point numbers or as small isolating intervals, or disks for complex roots.<br>
!>  Solving an equation \f$f(x) = g(x)\f$ is the same as finding the roots of the function \f$h(x) = f(x) – g(x)\f$.<br>
!>  Thus root-finding algorithms allow solving any equation defined by continuous functions.<br>
!>  However, most root-finding algorithms do not guarantee that they will find all the roots.<br>
!>  If such an algorithm does not find any root, it does not mean that no root exists.<br>
!>
!>  Most numerical root-finding methods use **iteration**, producing a sequence of numbers that hopefully converges towards the root as its limit.<br>
!>  They require one or more initial guesses of the root as starting values, then each iteration of the algorithm produces a successively more accurate approximation to the root.<br>
!>  Since the iteration must be stopped at some point, these methods produce an approximation to the root, not an exact solution.<br>
!>
!>  Bracketing methods
!>  ==================
!>
!>  Bracketing methods determine successively smaller intervals (brackets) that contain a root.<br>
!>  When the interval is small enough, then a root has been found.<br>
!>  They generally use the intermediate value theorem, which asserts that if a continuous function has values of
!>  opposite signs at the end points of an interval, then the function has at least one root in the interval.<br>
!>  Therefore, they *require to start with an interval such that the function takes opposite signs at the end points of the interval*.<br>
!>  There are [other methods](https://en.wikipedia.org/wiki/Descartes%27_rule_of_signs) for getting information on the number of roots of polynomials in an interval.<br>
!>
!>  Bisection method
!>  ----------------
!>
!>  The Bisection method consists of repeatedly bisecting a pre-specified interval known to contain at least one root.<br>
!>  The method selects the subinterval in which the function changes sign, and therefore must contain a root.<br>
!>  It is a very simple and robust, but relatively slow root-finding method.<br>
!>  As such, it is frequently used to obtain a rough approximation to the solution of a given root-finding problem.<br>
!>  The approximate solution is then used as a starting point for more rapidly converging methods.<br>
!>
!>  <b>The Bisection Algorithm</b><br>
!>  The Bisection method numerically solves the real-valued equation \f$f(x) = 0\f$.<br>
!>  The continuous function \f$f\f$ is defined on a search interval \f$[a, b]\f$ with \f$f(a)\f$ and \f$f(b)\f$ having opposite signs.<br>
!>  In such a case \f$a\f$ and \f$b\f$ are said to **bracket a root** \f$f\f$ must have at least one root in the search interval \f$(a, b)\f$.<br>
!>  At each step, the method divides the interval in two parts (halves) by computing the midpoint \f$c = (a + b) / 2\f$ of the interval and \f$f(c)\f$.<br>
!>  If \f$c\f$ is a root then the process has succeeded and stops.<br>
!>  Otherwise, there are now only two possibilities:<br>
!>  <ol>
!>      <li>    \f$f(a)\f$ and \f$f(c)\f$ have opposite signs and bracket a root.<br>
!>      <li>    \f$f(c)\f$ and \f$f(b)\f$ have opposite signs and bracket a root.<br>
!>  </ol>
!>  If the function has the same sign at the endpoints of an interval, the endpoints may or may not bracket roots of the function.<br>
!>  The Bisection method selects the subinterval that is guaranteed to be a bracket as the new interval to be used in the next step.<br>
!>  In this way an interval that contains a zero of \f$f\f$ is reduced in width by half at each step.<br>
!>  The process is continued until the interval is sufficiently small.<br>
!>
!>  \remark
!>  The finite precision of computers can cause serious problems for convergence of the Bisection method.<br>
!>  As such, the implementations of the method frequently include convergence tests or limits to the number of iterations.<br>
!>  Although \f$f\f$ is continuous, the finite precision of computer can preclude a function value ever being zero.<br>
!>  Additionally, the difference between the search interval limits cannot be less than the floating point precision of the computer.<br>
!>
!>  \note
!>  The Bisection method is also known as the **Interval-Halving method**, the **Binary Search method**, or the **Dichotomy method**.
!>
!>  False position (regula falsi) method
!>  ------------------------------------
!>
!>  The False Position method is a root-finding algorithm is a very old method for solving an equation with one unknown.<br>
!>  In simple terms, the method is the trial and error technique of using test (*false*) values for the variable and
!>  then adjusting the test value according to the outcome.<br>
!>  This is sometimes also referred to as *guess and check*.<br>
!>  Versions of the method predate the advent of algebra and the use of equations.<br>
!>
!>  <b>The False Position Algorithm</b><br>
!>
!>  The Regula Falsi method calculates the new solution estimate as the x-intercept of the line segment joining the endpoints of the function on the current bracketing interval of the root.<br>
!>  Essentially, the root is being approximated by replacing the actual function by a line segment on the bracketing interval and
!>  then using the classical double false position formula on that line segment.<br>
!>
!>  \image html pm_mathRoot@false.png width=500
!>
!>  Suppose that in the \f$k\f$-th iteration the bracketing interval is \f$(a_k, b_k)\f$ as illustrated above.<br>
!>  Construct the line through the points \f$(a_k, f(a_k))\f$ and \f$(b_k, f(b_k))\f$.<br>
!>  The equation of the line is given by,<br>
!>  \f{equation}{
!>      \large
!>      y - f(b_k) = \frac{f(b_k) - f(a_k)}{b_k - a_k} (x - b_k) ~.
!>  \f}
!>  Now choose \f$c_k\f$ to be the x-intercept of this line, that is, the value of \f$x\f$ for which \f$y = 0\f$, and substitute these values to obtain,<br>
!>  \f{equation}{
!>      \large
!>      f(b_k) + \frac{f(b_k) - f(a_k)}{b_k - a_k}(c_k - b_k) = 0 ~.
!>  \f}
!>  Solving this equation for \f$c_k\f$ yields,<br>
!>  \f{equation}{
!>      \large
!>      c_k = b_k - f(b_k) \frac{b_k - a_k}{f(b_k) - f(a_k)} = \frac{a_k f(b_k) - b_k f(a_k)}{f(b_k) - f(a_k)} ~.
!>  \f}
!>  The last symmetrical form has a computational advantage;<br>
!>  As a solution is approached, \f$a_k\f$ and \f$b_k\f$ will be very close together and nearly always of the same sign.<br>
!>  Such a subtraction can lose significant digits.<br>
!>  Because \f$f(b_k)\f$ and \f$f(a_k)\f$ are always of opposite sign, the subtraction in the numerator of the improved formula is effectively an addition (as is the subtraction in the denominator).<br>
!>  At iteration number \f$k\f$, the number \f$c_k\f$ is calculated as above and then, if \f$f(a_k)\f$ and \f$f(c_k)\f$ have the same sign, set \f$a_k + 1 = c_k\f$ and \f$b_k + 1 = b_k\f$.<br>
!>  Otherwise set \f$a_k + 1 = a_k\f$ and \f$b_k + 1 = c_k\f$.<br>
!>  This process is repeated until the root is approximated sufficiently well.<br>
!>  The above formula is also used in the Secant method.<br>
!>  However, the Secant method always retains the last two computed points.<br>
!>  Therefore, while the Secant method is slightly faster, it does not preserve the root bracketing and may not converge.<br>
!>  The fact that **the Regula Falsi root-finding method always converges** makes it a good choice when speed is needed.<br>
!>  However, its rate of convergence can drop below that of the bisection method.<br>
!>
!>  <b>The Convergence of the False Position Algorithm</b><br>
!>
!>  The convergence rate of the False Position method can be better than the Bisection method but is worse than the Secant method.<br>
!>  The method often has a superlinear convergence rate.<br>
!>  However, estimation of the exact order of convergence is difficult (e.g., see Numerical Recipes in Fortran by Press et al. 1992 for a discussion).<br>
!>
!>  \remark
!>  The finite precision of computers can cause problems for convergence of the False Position method.<br>
!>  As such, the implementations of the method frequently include convergence tests or limits to the number of iterations.<br>
!>  Although \f$f\f$ is continuous, the finite precision of computer can preclude a function value ever being zero.<br>
!>  Additionally, the difference between the search interval limits cannot be less than the floating point precision of the computer.<br>
!>
!>  Iterative methods
!>  =================
!>
!>  The Secant method
!>  -----------------
!>
!>  The Secant method is a root-finding algorithm that uses a succession of roots of secant lines to better approximate a root of a function.<br>
!>  The Secant method can be thought of as a finite-difference approximation of the [Newton method](@ref pm_mathRoot).<br>
!>  However, the Secant method predates the Newton method by over 3000 years.<br>
!>
!>  <b>The Secant Algorithm</b><br>
!>  The Secant method is defined by the following recursive relation,<br>
!>  \f{eqnarray}{
!>      \large
!>      x_{n}
!>      &=& x_{n-1} - f(x_{n-1}) \frac{x_{n-1} - x_{n-2}} {f(x_{n-1}) - f(x_{n-2})} ~, \\
!>      &=& \frac {x_{n-2} f(x_{n-1}) - x_{n-1} f(x_{n-2})} {f(x_{n-1}) - f(x_{n-2})} ~,
!>  \f}
!>  where the two initial values \f$x_0\f$ and \f$x_1\f$ should be chosen close to the desired zero.<br>
!>
!>  \remark
!>  Unlike the other iterative root-finding methods such as the [Bisection](@ref pm_mathRoot) or the [Brent](@ref pm_mathRoot) methods,
!>  the Secant method does **not** require the initial starting values \f$x_0\f$ and \f$x_1\f$ to bracket the root of the function.
!>
!>  <b>Convergence of the Secant Algorithm</b><br>
!>  The iterates \f$x_n\f$ of the Secant method converge to a root of a continuous function if the initial values \f$x_0\f$ and \f$x_1\f$ are sufficiently close to the root.<br>
!>  The order of convergence is the Golden Ratio \f$\phi = 1.618\f$, so that the limiting error in the root is,<br>
!>  \f{equation}{
!>      \large
!>      \lim_{k\rightarrow+\infty} |\epsilon_{k}| \propto |\epsilon_{k - 1}|^{1.618} ~.
!>  \f}
!>  The rate of convergence is therefore faster than the [Bisection](@ref pm_mathRoot).<br>
!>  In particular, the convergence is **super-linear**, but **sub-quadratic**.<br>
!>  This convergence rate only holds under some technical conditions, namely the **function must be twice continuously differentiable and the root in question be simple** (i.e., with multiplicity 1).<br>
!>  If the initial values are not close enough to the root, then there is no guarantee that the Secant method converges.<br>
!>  There is no general definition of *close enough*, but the criterion has to do with how *wiggly* the function is on the interval \f$[x_0, x_1]\f$.<br>
!>  For example, if the function is differentiable on the interval and there is a point where \f$f′(x) = 0\f$ on the interval, then the algorithm may not converge.<br>
!>  For functions that are not sufficiently continuous, the algorithm can therefore not be guaranteed to converge: Local behavior might send it off towards infinity.<br>
!>
!>  \remark
!>  The finite precision of computers can cause problems for convergence of the Secant method.<br>
!>  As such, the implementations of the method frequently include convergence tests or limits to the number of iterations.<br>
!>  Although \f$f\f$ is continuous, the finite precision of computer can preclude a function value ever being zero.<br>
!>  Additionally, the difference between the search interval limits cannot be less than the floating point precision of the computer.<br>
!>
!>  The Newton method
!>  -----------------
!>
!>  The Newton method, also known as the Newton–Raphson method, named after Isaac Newton and Joseph Raphson,
!>  is a root-finding algorithm which produces successively better approximations to the roots (or zeroes) of a real-valued function.<br>
!>
!>  <b>The Newton Algorithm</b><br>
!>
!>  The most basic version of the method starts with,<br>
!>  <ol>
!>      <li>    a single-variable function \f$f(x)\f$ defined for a real variable \f$x\f$,<br>
!>      <li>    the function derivative \f$f′(x)\f$,<br>
!>      <li>    an initial guess \f$x_0\f$ for a root of \f$f(x)\f$.<br>
!>  </ol>
!>  If the function satisfies sufficient assumptions and the initial guess is close, then,<br>
!>  \f{equation}{
!>      \large
!>      x_{1} = x_{0} - \frac{f(x_{0})}{f'(x_{0})} ~,
!>  \f}
!>  is a better approximation of the root than \f$x_0\f$ (as illustrated below).<br>
!>
!>  \image html pm_mathRoot@newton.gif width=500
!>
!>  Geometrically, \f$(x_1, 0)\f$ is the intersection of the x-axis and the tangent of the graph of \f$f\f$ at \f$(x_0, f(x_0))\f$.<br>
!>  In other words, the improved guess is the unique root of the linear approximation at the initial point.<br>
!>  The process is repeated as,<br>
!>  \f{equation}{
!>      \large
!>      x_{n+1} = x_{n} - \frac {f(x_{n})}{f'(x_{n})} ~,
!>  \f}
!>  until a sufficiently precise value is reached.<br>
!>  The Newton method can also be extended to complex functions and to systems of equations.<br>
!>
!>  <b>The Convergence of the Newton Algorithm</b><br>
!>
!>  The Newton method will usually converge, provided the initial guess is *close enough* to the unknown root of the function and \f$f\prime(x_0) \neq 0\f$.<br>
!>  Furthermore, for a root of multiplicity \f$1\f$, the convergence is at least quadratic in a neighborhood of the root.<br>
!>  This intuitively means that the number of correct digits roughly doubles in every step.<br>
!>  More details can be found in the analysis section below.
!>
!>  <b>Practical Considerations for the Newton Algorithm</b><br>
!>
!>  The Newton method is a powerful technique.<br>
!>  Convergence is generally quadratic.<br>
!>  However, there are some difficulties associated with the method.<br>
!>  <ul>
!>      <li>    **The Newton method requires the function derivative.**<br>
!>              The Newton method requires that the derivative can be calculated directly.<br>
!>              An analytical expression for the derivative may not be easily obtainable or could be expensive to evaluate.<br>
!>              In these situations, it may be appropriate to approximate the derivative by using the slope of a line through two nearby points on the function.<br>
!>              However, numerical approximation of the derivative degrade the convergence rate and performance to the Secant method which is slower than the Newton method.<br>
!>      <li>    **The Newton method can fail to converge to the root.**<br>
!>              The proof quadratic convergence of the Newton method requires certain assumptions that may not always hold.<br>
!>              One should therefore review the assumptions before using the Newton method.<br>
!>              For situations where the method fails to converge, it is because the assumptions made in the proof convergence are not met.<br>
!>      <li>    **The Newton method can overshoot.**<br>
!>              If the first derivative is not well behaved in the neighborhood of a particular root, the method may overshoot, and diverge from that root.<br>
!>              An example of a function with one root, for which the derivative is not well behaved in the neighborhood of the root, is
!>              \f{equation}{
!>                  \large
!>                  f(x) = |x|^{a} ~, \quad 0 < a < \frac{1}{2} ~,
!>              \f}
!>              for which the root will be overshot and the sequence of \f$x\f$ will diverge.<br>
!>              For \f$a = 1/2\f$, the root will still be overshot, but the sequence will oscillate between two values.<br>
!>              For \f$1/2 < a < 1\f$, the root will still be overshot but the sequence will converge.<br>
!>              For \f$a \geq 1\f$ the root will not be overshot at all.<br>
!>              In some cases, the Newton method can be stabilized by using successive over-relaxation.<br>
!>      <li>    **The Newton method can prematurely terminate upon encountering stationary points of the function.**<br>
!>              If a stationary point of the function is encountered, the derivative is zero and the method will terminate due to division by zero.<br>
!>      <li>    **The Newton method can fail due to poor initial starting point for the search.**<br>
!>              A large error in the initial estimate can contribute to non-convergence of the algorithm.<br>
!>              To overcome this problem one can linearize the function that is being optimized.<br>
!>              Good initial estimates lie close to the final globally optimal parameter estimate.<br>
!>  </ul>
!>
!>  \remark
!>  <b>The procedures of this module use a hybrid Newton-Bisection method to resolve the difficulties mentioned above.</b><br>
!>  See [this article](https://en.wikipedia.org/wiki/Newton%27s_method#Failure_analysis) for an exhaustive discussion of the failures and remedies.<br>
!>
!>  \remark
!>  The finite precision of computers can cause problems for convergence of the Newton method.<br>
!>  As such, the implementations of the method frequently include convergence tests or limits to the number of iterations.<br>
!>  Although \f$f\f$ is continuous, the finite precision of computer can preclude a function value ever being zero.<br>
!>  Additionally, the difference between the search interval limits cannot be less than the floating point precision of the computer.<br>
!>
!>  The Halley method
!>  -----------------
!>
!>  In numerical analysis, the Halley method is a root-finding algorithm used for functions of one real variable **with a continuous second derivative**.<br>
!>  It is named after its inventor [Edmond Halley](https://en.wikipedia.org/wiki/Edmond_Halley).<br>
!>  The algorithm is second in the class of Householder methods, after Newton method.<br>
!>  Like the latter, it iteratively produces a sequence of approximations to the root.<br>
!>  The **rate of convergence** to the root is **cubic**.<br>
!>  The Halley method exactly finds the roots of a linear-over-linear Padé approximation to the function,
!>  in contrast to the Newton method or the Secant method which approximate the function linearly,
!>  or the Muller method which approximates the function quadratically.<br>
!>
!>  <b>The Halley Algorithm</b><br>
!>
!>  The Halley method solves the nonlinear equation \f$f(x) = 0\f$.<br>
!>  In this case, the function \f$f\f$ has to be a function of one real variable.<br>
!>  The method consists of a sequence of iterations:<br>
!>  \f{equation}{
!>      x_{{n+1}} = x_{n} - \frac{2f(x_{n})f'(x_{n})}{2{[f'(x_{n})]}^{2}-f(x_{n})f''(x_{n})} ~,
!>  \f}
!>  beginning with an initial guess \f$x_0\f$.<br>
!>
!>  <b>The Convergence of the Halley Algorithm</b><br>
!>
!>  The Halley method converges **cubically** to the function root.<br>
!>  That is, each iteration triples the number of significant digits in the final root.<br>
!>  In comparison, two steps of Newton-Raphson quadruple the number of digits in the final root.<br>
!>
!>  <b>Practical Considerations for the Halley Algorithm</b><br>
!>
!>  The use of the Halley method is sensible only when it is easy to calculate the second derivative of the target function.<br>
!>  The basin of convergence of the Halley method is not guaranteed to be larger than that of the Newton method.<br>
!>  This means that the Halley method does not necessarily yield faster convergence rates as evidenced by the examples of the generic interfaces of this module.<br>
!>  he current implementation of the Halley method in this module resolves over-compensations
!>  by the second derivative (leading to wrong search directions) by reverting the search method to a Newton-Raphson step.<br>
!>  If the method sends the search out of the initial user-specified bracket of the root, then the algorithm falls back to the bisection method to correct the search.<br>
!>
!>  The Schroder method
!>  -------------------
!>
!>  Similar to the Halley method, the Schroder method solves the nonlinear equation \f$f(x) = 0\f$ using the second derivative.<br>
!>  However, unlike the Newton and the Halley methods, the Schroder method is known to work well in the presence of multiple roots.<br>
!>  The method consists of a sequence of iterations:<br>
!>  \f{equation}{
!>      x_{{n+1}} = x_{n} - \frac{f(x_n)}{f'(x_n)} - \frac{[f(x_n)]^2 f''(x_n)}{2[f'(x_n)]^3} ~,
!>  \f}
!>  beginning with an initial guess \f$x_0\f$.<br>
!>  See [Stewart, 1993, G. W. On Infinitely Many Algorithms for Solving Equations](http://drum.lib.umd.edu/handle/1903/577)
!>  for the English translation of the original paper of Schroder and the derivation of the above equation on page 13.<br>
!>
!>  Hybrid methods
!>  ==============
!>
!>  The Brent method
!>  ----------------
!>
!>  The Brent method is a **hybrid** root-finding algorithm combining<br>
!>  <ul>
!>      <li>    the **Bisection method**,<br>
!>      <li>    the **Secant method**, and<br>
!>      <li>    the **inverse quadratic interpolation**.<br>
!>  </ul>
!>  It has the **reliability of the Bisection method** but potentially as fast as some of the less-reliable methods above.<br>
!>  The algorithm tries to use the potentially **fast-converging** Secant method or inverse quadratic interpolation if possible.<br>
!>  If necessary, it falls back to the more robust Bisection method.<br>
!>  The Brent method is due to [Richard Brent](https://en.wikipedia.org/wiki/Richard_P._Brent) and builds on an earlier algorithm by [Theodorus Dekker](https://en.wikipedia.org/wiki/Theodorus_Dekker).<br>
!>
!>  The Quadratic method employed in the Brent method works well only when the function behaves smoothly.<br>
!>  However, they run the serious risk of giving very bad estimates of the root or causing machine failure by an inappropriate division by a very small number.<br>
!>  The Brent method guards against this problem by maintaining the search brackets on the root and checking where the interpolation would land before carrying out divisions.<br>
!>  When the corrections due to the Quadratic method would not land within the search bounds, or when the bounds are not collapsing rapidly enough, the algorithm takes a bisection step.<br>
!>  Thus, **the Brent method combines the sureness of the Bisection method with the speed of a higher-order method when appropriate**.<br>
!>
!>  \remark
!>  The brent algorithm implemented in this module is a reimplementation of the original [FORTRAN77 Brent method from the NetLib library](https://netlib.org/go/zeroin.f)
!>  combined with improvements inspired by the [Numerical Recipes in Fortran, Press et al. 1991](http://numerical.recipes/).<br>
!>
!>  \note
!>  The algorithm is also known as the **Van Wijngaarden–Dekker–Brent root-finding method** or the **Brent–Dekker root-finding method**.<br>
!>
!>  \remark
!>  The finite precision of computers can cause problems for convergence of the Brent method.<br>
!>  As such, the implementations of the method frequently include convergence tests or limits to the number of iterations.<br>
!>  Although \f$f\f$ is continuous, the finite precision of computer can preclude a function value ever being zero.<br>
!>  Additionally, the difference between the search interval limits cannot be less than the floating point precision of the computer.<br>
!>
!>  The Ridders method
!>  ------------------
!>
!>  The Ridders method is a root-finding algorithm based on the False-Position method
!>  and the use of an exponential function to successively approximate a root of a continuous function.<br>
!>  The method is due to C. Ridders.<br>
!>
!>  <b>The Ridders Algorithm</b><br>
!>
!>  Given two values of the independent variable, \f$x_0 < \mathrm{root} < x_2\f$, the method begins by evaluating the function at the midpoint \f$x_1 = (x_0 + x_2) / 2\f$.<br>
!>  One then finds the unique exponential function \f$e^{ax}\f$ such that the function \f$h(x) = f(x) \exp(ax)\f$ satisfies \f$h(x_1) = (h(x_0) + h(x_2)) / 2\f$.<br>
!>  Specifically, the parameter \f$a\f$ is determined by,<br>
!>  \f{equation}{
!>      \large
!>      \exp\big(a(x_1 - x_0)\big) = \frac{f(x_1) - \mathrm{sign}\big(f(x_0)\big){\sqrt{f(x_1)^2 - f(x_0) f(x_2)}}}{f(x_2)} ~.
!>  \f}
!>  The False-Position method is then applied to the points \f$(x_0, h(x_0))\f$ and \f$(x_2, h(x_2))\f$,
!>  leading to a new value \f$x_3\f$ between \f$x_0\f$ and \f$x_2\f$,
!>  \f{equation}{
!>      \large
!>      x_3 = x_1 + (x_1 - x_0) \frac{\mathrm{sign}\big(f(x_0)\big) f(x_1)} {\sqrt{f(x_1)^2 - f(x_0) f(x_2)}} ~,
!>  \f}
!>  which will be used as one of the two bracketing values in the next step of the iteration.<br>
!>  The other bracketing value is taken to be \f$x_1\f$ if \f$f(x_1) f(x_3) < 0\f$ (well-behaved case),
!>  or otherwise whichever of \f$x_0\f$ and \f$x_2\f$ has function value of opposite sign to \f$f(x_3)\f$.<br>
!>  The procedure terminates when the desired accuracy is achieved.<br>
!>
!>  <b>The Convergence of the Ridders Algorithm</b><br>
!>
!>  The Ridders method is simpler than the Muller method or the Brent method but with similar performance.<br>
!>  The algorithm converges quadratically when the function is well-behaved.<br>
!>  This implies that the number of additional significant digits found at each step approximately doubles.<br>
!>  However, since the function has to be evaluated twice per iteration, the overall order of convergence of the method is \f$\sqrt{2}\f$.<br>
!>  The root remains bracketed and the length of the bracketing interval at least halves on each iteration, even the function is not well-behaved.<br>
!>  In all circumstances, the convergence is guaranteed.
!>
!>  \remark
!>  The finite precision of computers can cause problems for convergence of the Ridders method.<br>
!>  As such, the implementations of the method frequently include convergence tests or limits to the number of iterations.<br>
!>  Although \f$f\f$ is continuous, the finite precision of computer can preclude a function value ever being zero.<br>
!>  Additionally, the difference between the search interval limits cannot be less than the floating point precision of the computer.<br>
!>
!>  The TOMS748 method
!>  ------------------
!>
!>  The [TOMS748 algorithm](https://na.math.kit.edu/alefeld/download/1995_Algorithm_748_Enclosing_Zeros_of_Continuous_Functions.pdf)
!>  is a hybrid root-finding algorithm introduced in,<br>
!>  <ul>
!>      <li>    Alefeld, G. E., Potra, F. A., Shi, Yixun (1995). "Algorithm 748: Enclosing Zeros of Continuous Functions".
!>              ACM Transactions on Mathematical Software. 21 (3): 327–344. doi:10.1145/210089.210111
!>  </ul>
!>  that uses a mixture of cubic, quadratic and linear (secant) interpolation to locate the root of a function.<br>
!>  While there is widespread claims of the superior performance of this algorithm compared to the Brent method,
!>  there are [counter arguments](https://www.nongnu.org/lmi/toms_748.html) for such benchmarks.<br>
!>
!>  The optimal root-finding method
!>  ===============================
!>
!>  The best root-finding algorithm depends on the problem being solved:<br>
!>  <ol>
!>      <li>    The Brent method is the recommended **algorithm of choice**
!>              for general one-dimensional root-finding problems **without** knowledge of function derivative (gradient).<br>
!>      <li>    The [Newton-Raphson method](@ref pm_mathRoot) is the recommended **algorithm of choice**
!>              for general one-dimensional root-finding problems **with** knowledge of function derivative (gradient).<br>
!>  </ol>
!>
!>  \see
!>  [pm_arraySearch](@ref pm_arraySearch)<br>
!>  [getPolyRoot](@ref pm_polynomial::getPolyRoot)<br>
!>  [setPolyRoot](@ref pm_polynomial::setPolyRoot)<br>
!>  [Regula Falsi](https://en.wikipedia.org/wiki/Regula_falsi)<br>
!>  [The Secant Algorithm](https://en.wikipedia.org/wiki/Secant_method)<br>
!>  [The Bisection Algorithm](https://en.wikipedia.org/wiki/Bisection_method)<br>
!>  [Root-Finding Algorithms](https://en.wikipedia.org/wiki/Root-finding_algorithms)<br>
!>  [Root-Finding Algorithms](https://en.wikipedia.org/wiki/Root-finding_algorithms)<br>
!>  [Method of False Position](https://mathworld.wolfram.com/MethodofFalsePosition.html#:~:text=An%20algorithm%20for%20finding%20roots,1992)<br>
!>  Brent, Richard P., 1971, *An algorithm with guaranteed convergence for finding a zero of a function*, The computer journal, 14, 4, 422-425<br>
!>  Newton, Richard P., 1971, *An algorithm with guaranteed convergence for finding a zero of a function*, The computer journal, 14, 4, 422-425<br>
!>  Ridders, C. F. J. *A New Algorithm for Computing a Single Root of a Real Continuous Function.* IEEE Trans. Circuits Systems 26, 979-980, 1979<br>
!>  Schröder, E. "Über unendlich viele Algorithmen zur Auflösung der Gleichungen." Math. Ann. 2, 317-365, 1870<br>
!>  [RootsFortran](https://github.com/jacobwilliams/roots-fortran): An extensive Fortran library for root-finding by [Jacob Williams](https://github.com/jacobwilliams)<br>
!>  Numerical Recipes in Fortran, Press et al. 1992<br>
!>  [SLATEC cdzro.f](https://netlib.org/slatec/src/)<br>
!>  [Boost library](https://www.boost.org/doc/libs/1_82_0/libs/math/doc/html/root_finding.html)<br>
!>
!>  \test
!>  [test_pm_mathRoot](@ref test_pm_mathRoot)
!>
!>  \todo
!>  \pmed
!>  The [Muller method](https://en.wikipedia.org/wiki/Muller%27s_method) of root-finding should be implemented.<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_mathRoot

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathRoot"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is an `abstract` derived type for constructing concrete derived types to
    !>  distinguish various procedure signatures that require root-finding methods (e.g., Bisection, False Position, Secant, Newton, Brent, Ridders, ...).<br>
    !>
    !>  \details
    !>  This `abstract` derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users must use `parameter` objects instantiated from the concrete subclasses of this parent `abstract` derived type.<br>
    !>
    !>  \see
    !>  [brent](@ref pm_mathRoot::brent)<br>
    !>  [false](@ref pm_mathRoot::false)<br>
    !>  [secant](@ref pm_mathRoot::secant)<br>
    !>  [halley](@ref pm_mathRoot::halley)<br>
    !>  [newton](@ref pm_mathRoot::newton)<br>
    !>  [ridders](@ref pm_mathRoot::ridders)<br>
    !>  [toms748](@ref pm_mathRoot::toms748)<br>
    !>  [schroder](@ref pm_mathRoot::schroder)<br>
    !>  [bisection](@ref pm_mathRoot::bisection)<br>
    !>  [brent_type](@ref pm_mathRoot::brent_type)<br>
    !>  [false_type](@ref pm_mathRoot::false_type)<br>
    !>  [secant_type](@ref pm_mathRoot::secant_type)<br>
    !>  [halley_type](@ref pm_mathRoot::halley_type)<br>
    !>  [newton_type](@ref pm_mathRoot::newton_type)<br>
    !>  [ridders_type](@ref pm_mathRoot::ridders_type)<br>
    !>  [toms748_type](@ref pm_mathRoot::toms748_type)<br>
    !>  [schroder_type](@ref pm_mathRoot::schroder_type)<br>
    !>  [bisection_type](@ref pm_mathRoot::bisection_type)<br>
    !>  [iteration_type](@ref pm_mathRoot::iteration_type)<br>
    !>  [bracket_type](@ref pm_mathRoot::bracket_type)<br>
    !>  [hybrid_type](@ref pm_mathRoot::hybrid_type)<br>
    !>  [method_type](@ref pm_mathRoot::method_type)<br>
    !>
    !>  \final{method_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, abstract :: method_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is an `abstract` derived type for constructing concrete derived types to
    !>  distinguish various procedure signatures that require **bracketing** root-finding methods (e.g., Bisection, False Position, ...).<br>
    !>
    !>  \details
    !>  This `abstract` derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users must use `parameter` objects instantiated from the concrete subclasses of this parent `abstract` derived type.<br>
    !>
    !>  \see
    !>  [brent](@ref pm_mathRoot::brent)<br>
    !>  [false](@ref pm_mathRoot::false)<br>
    !>  [secant](@ref pm_mathRoot::secant)<br>
    !>  [halley](@ref pm_mathRoot::halley)<br>
    !>  [newton](@ref pm_mathRoot::newton)<br>
    !>  [ridders](@ref pm_mathRoot::ridders)<br>
    !>  [toms748](@ref pm_mathRoot::toms748)<br>
    !>  [schroder](@ref pm_mathRoot::schroder)<br>
    !>  [bisection](@ref pm_mathRoot::bisection)<br>
    !>  [brent_type](@ref pm_mathRoot::brent_type)<br>
    !>  [false_type](@ref pm_mathRoot::false_type)<br>
    !>  [secant_type](@ref pm_mathRoot::secant_type)<br>
    !>  [halley_type](@ref pm_mathRoot::halley_type)<br>
    !>  [newton_type](@ref pm_mathRoot::newton_type)<br>
    !>  [ridders_type](@ref pm_mathRoot::ridders_type)<br>
    !>  [toms748_type](@ref pm_mathRoot::toms748_type)<br>
    !>  [schroder_type](@ref pm_mathRoot::schroder_type)<br>
    !>  [bisection_type](@ref pm_mathRoot::bisection_type)<br>
    !>  [iteration_type](@ref pm_mathRoot::iteration_type)<br>
    !>  [bracket_type](@ref pm_mathRoot::bracket_type)<br>
    !>  [hybrid_type](@ref pm_mathRoot::hybrid_type)<br>
    !>  [method_type](@ref pm_mathRoot::method_type)<br>
    !>
    !>  \final{bracket_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, abstract, extends(method_type) :: bracket_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is an `abstract` derived type for constructing concrete derived types to
    !>  distinguish various procedure signatures that require **iterative** root-finding methods (e.g., Secant, Newton, ...).<br>
    !>
    !>  \details
    !>  This `abstract` derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users must use `parameter` objects instantiated from the concrete subclasses of this parent `abstract` derived type.<br>
    !>
    !>  \see
    !>  [brent](@ref pm_mathRoot::brent)<br>
    !>  [false](@ref pm_mathRoot::false)<br>
    !>  [secant](@ref pm_mathRoot::secant)<br>
    !>  [halley](@ref pm_mathRoot::halley)<br>
    !>  [newton](@ref pm_mathRoot::newton)<br>
    !>  [ridders](@ref pm_mathRoot::ridders)<br>
    !>  [toms748](@ref pm_mathRoot::toms748)<br>
    !>  [schroder](@ref pm_mathRoot::schroder)<br>
    !>  [bisection](@ref pm_mathRoot::bisection)<br>
    !>  [brent_type](@ref pm_mathRoot::brent_type)<br>
    !>  [false_type](@ref pm_mathRoot::false_type)<br>
    !>  [secant_type](@ref pm_mathRoot::secant_type)<br>
    !>  [halley_type](@ref pm_mathRoot::halley_type)<br>
    !>  [newton_type](@ref pm_mathRoot::newton_type)<br>
    !>  [ridders_type](@ref pm_mathRoot::ridders_type)<br>
    !>  [toms748_type](@ref pm_mathRoot::toms748_type)<br>
    !>  [schroder_type](@ref pm_mathRoot::schroder_type)<br>
    !>  [bisection_type](@ref pm_mathRoot::bisection_type)<br>
    !>  [iteration_type](@ref pm_mathRoot::iteration_type)<br>
    !>  [bracket_type](@ref pm_mathRoot::bracket_type)<br>
    !>  [hybrid_type](@ref pm_mathRoot::hybrid_type)<br>
    !>  [method_type](@ref pm_mathRoot::method_type)<br>
    !>
    !>  \final{iteration_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, abstract, extends(method_type) :: iteration_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is an `abstract` derived type for constructing concrete derived types to
    !>  distinguish various procedure signatures that require **iterative** root-finding methods (e.g., Secant, Newton, ...).<br>
    !>
    !>  \details
    !>  This `abstract` derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users must use `parameter` objects instantiated from the concrete subclasses of this parent `abstract` derived type.<br>
    !>
    !>  \see
    !>  [brent](@ref pm_mathRoot::brent)<br>
    !>  [false](@ref pm_mathRoot::false)<br>
    !>  [secant](@ref pm_mathRoot::secant)<br>
    !>  [halley](@ref pm_mathRoot::halley)<br>
    !>  [newton](@ref pm_mathRoot::newton)<br>
    !>  [ridders](@ref pm_mathRoot::ridders)<br>
    !>  [toms748](@ref pm_mathRoot::toms748)<br>
    !>  [schroder](@ref pm_mathRoot::schroder)<br>
    !>  [bisection](@ref pm_mathRoot::bisection)<br>
    !>  [brent_type](@ref pm_mathRoot::brent_type)<br>
    !>  [false_type](@ref pm_mathRoot::false_type)<br>
    !>  [secant_type](@ref pm_mathRoot::secant_type)<br>
    !>  [halley_type](@ref pm_mathRoot::halley_type)<br>
    !>  [newton_type](@ref pm_mathRoot::newton_type)<br>
    !>  [ridders_type](@ref pm_mathRoot::ridders_type)<br>
    !>  [toms748_type](@ref pm_mathRoot::toms748_type)<br>
    !>  [schroder_type](@ref pm_mathRoot::schroder_type)<br>
    !>  [bisection_type](@ref pm_mathRoot::bisection_type)<br>
    !>  [iteration_type](@ref pm_mathRoot::iteration_type)<br>
    !>  [bracket_type](@ref pm_mathRoot::bracket_type)<br>
    !>  [hybrid_type](@ref pm_mathRoot::hybrid_type)<br>
    !>  [method_type](@ref pm_mathRoot::method_type)<br>
    !>
    !>  \final{hybrid_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, abstract, extends(method_type) :: hybrid_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the use of the Brent method of root-finding.<br>
    !>
    !>  \details
    !>  See the root-finding section in the documentation of [pm_mathRoot](@ref pm_mathRoot) for more details about this root-finding method.<br>
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [brent](@ref pm_mathRoot::brent) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [brent](@ref pm_mathRoot::brent)<br>
    !>  [false](@ref pm_mathRoot::false)<br>
    !>  [secant](@ref pm_mathRoot::secant)<br>
    !>  [halley](@ref pm_mathRoot::halley)<br>
    !>  [newton](@ref pm_mathRoot::newton)<br>
    !>  [ridders](@ref pm_mathRoot::ridders)<br>
    !>  [toms748](@ref pm_mathRoot::toms748)<br>
    !>  [schroder](@ref pm_mathRoot::schroder)<br>
    !>  [bisection](@ref pm_mathRoot::bisection)<br>
    !>  [brent_type](@ref pm_mathRoot::brent_type)<br>
    !>  [false_type](@ref pm_mathRoot::false_type)<br>
    !>  [secant_type](@ref pm_mathRoot::secant_type)<br>
    !>  [halley_type](@ref pm_mathRoot::halley_type)<br>
    !>  [newton_type](@ref pm_mathRoot::newton_type)<br>
    !>  [ridders_type](@ref pm_mathRoot::ridders_type)<br>
    !>  [toms748_type](@ref pm_mathRoot::toms748_type)<br>
    !>  [schroder_type](@ref pm_mathRoot::schroder_type)<br>
    !>  [bisection_type](@ref pm_mathRoot::bisection_type)<br>
    !>  [iteration_type](@ref pm_mathRoot::iteration_type)<br>
    !>  [bracket_type](@ref pm_mathRoot::bracket_type)<br>
    !>  [hybrid_type](@ref pm_mathRoot::hybrid_type)<br>
    !>  [method_type](@ref pm_mathRoot::method_type)<br>
    !>
    !>  \final{brent_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(hybrid_type) :: brent_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [brent_type](@ref pm_mathRoot::brent_type) that is exclusively used
    !>  to signify the use of Brent method of root-finding within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  See the root-finding section in the documentation of [pm_mathRoot](@ref pm_mathRoot) for more details about this root-finding method.<br>
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [brent](@ref pm_mathRoot::brent)<br>
    !>  [false](@ref pm_mathRoot::false)<br>
    !>  [secant](@ref pm_mathRoot::secant)<br>
    !>  [halley](@ref pm_mathRoot::halley)<br>
    !>  [newton](@ref pm_mathRoot::newton)<br>
    !>  [ridders](@ref pm_mathRoot::ridders)<br>
    !>  [toms748](@ref pm_mathRoot::toms748)<br>
    !>  [schroder](@ref pm_mathRoot::schroder)<br>
    !>  [bisection](@ref pm_mathRoot::bisection)<br>
    !>  [brent_type](@ref pm_mathRoot::brent_type)<br>
    !>  [false_type](@ref pm_mathRoot::false_type)<br>
    !>  [secant_type](@ref pm_mathRoot::secant_type)<br>
    !>  [halley_type](@ref pm_mathRoot::halley_type)<br>
    !>  [newton_type](@ref pm_mathRoot::newton_type)<br>
    !>  [ridders_type](@ref pm_mathRoot::ridders_type)<br>
    !>  [toms748_type](@ref pm_mathRoot::toms748_type)<br>
    !>  [schroder_type](@ref pm_mathRoot::schroder_type)<br>
    !>  [bisection_type](@ref pm_mathRoot::bisection_type)<br>
    !>  [iteration_type](@ref pm_mathRoot::iteration_type)<br>
    !>  [bracket_type](@ref pm_mathRoot::bracket_type)<br>
    !>  [hybrid_type](@ref pm_mathRoot::hybrid_type)<br>
    !>  [method_type](@ref pm_mathRoot::method_type)<br>
    !>
    !>  \final{brent}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(brent_type), parameter :: brent = brent_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: brent
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the use of the TOMS748 method of root-finding.<br>
    !>
    !>  \details
    !>  See the root-finding section in the documentation of [pm_mathRoot](@ref pm_mathRoot) for more details about this root-finding method.<br>
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [toms748](@ref pm_mathRoot::toms748) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [brent](@ref pm_mathRoot::brent)<br>
    !>  [false](@ref pm_mathRoot::false)<br>
    !>  [secant](@ref pm_mathRoot::secant)<br>
    !>  [halley](@ref pm_mathRoot::halley)<br>
    !>  [newton](@ref pm_mathRoot::newton)<br>
    !>  [ridders](@ref pm_mathRoot::ridders)<br>
    !>  [toms748](@ref pm_mathRoot::toms748)<br>
    !>  [schroder](@ref pm_mathRoot::schroder)<br>
    !>  [bisection](@ref pm_mathRoot::bisection)<br>
    !>  [brent_type](@ref pm_mathRoot::brent_type)<br>
    !>  [false_type](@ref pm_mathRoot::false_type)<br>
    !>  [secant_type](@ref pm_mathRoot::secant_type)<br>
    !>  [halley_type](@ref pm_mathRoot::halley_type)<br>
    !>  [newton_type](@ref pm_mathRoot::newton_type)<br>
    !>  [ridders_type](@ref pm_mathRoot::ridders_type)<br>
    !>  [toms748_type](@ref pm_mathRoot::toms748_type)<br>
    !>  [schroder_type](@ref pm_mathRoot::schroder_type)<br>
    !>  [bisection_type](@ref pm_mathRoot::bisection_type)<br>
    !>  [iteration_type](@ref pm_mathRoot::iteration_type)<br>
    !>  [bracket_type](@ref pm_mathRoot::bracket_type)<br>
    !>  [hybrid_type](@ref pm_mathRoot::hybrid_type)<br>
    !>  [method_type](@ref pm_mathRoot::method_type)<br>
    !>
    !>  \final{toms748_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(hybrid_type) :: toms748_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [toms748_type](@ref pm_mathRoot::toms748_type) that is exclusively used
    !>  to signify the use of TOMS748 method of root-finding within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  See the root-finding section in the documentation of [pm_mathRoot](@ref pm_mathRoot) for more details about this root-finding method.<br>
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [brent](@ref pm_mathRoot::brent)<br>
    !>  [false](@ref pm_mathRoot::false)<br>
    !>  [secant](@ref pm_mathRoot::secant)<br>
    !>  [halley](@ref pm_mathRoot::halley)<br>
    !>  [newton](@ref pm_mathRoot::newton)<br>
    !>  [ridders](@ref pm_mathRoot::ridders)<br>
    !>  [toms748](@ref pm_mathRoot::toms748)<br>
    !>  [schroder](@ref pm_mathRoot::schroder)<br>
    !>  [bisection](@ref pm_mathRoot::bisection)<br>
    !>  [brent_type](@ref pm_mathRoot::brent_type)<br>
    !>  [false_type](@ref pm_mathRoot::false_type)<br>
    !>  [secant_type](@ref pm_mathRoot::secant_type)<br>
    !>  [halley_type](@ref pm_mathRoot::halley_type)<br>
    !>  [newton_type](@ref pm_mathRoot::newton_type)<br>
    !>  [ridders_type](@ref pm_mathRoot::ridders_type)<br>
    !>  [toms748_type](@ref pm_mathRoot::toms748_type)<br>
    !>  [schroder_type](@ref pm_mathRoot::schroder_type)<br>
    !>  [bisection_type](@ref pm_mathRoot::bisection_type)<br>
    !>  [iteration_type](@ref pm_mathRoot::iteration_type)<br>
    !>  [bracket_type](@ref pm_mathRoot::bracket_type)<br>
    !>  [hybrid_type](@ref pm_mathRoot::hybrid_type)<br>
    !>  [method_type](@ref pm_mathRoot::method_type)<br>
    !>
    !>  \final{toms748}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(toms748_type), parameter :: toms748 = toms748_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: toms748
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the use of the False-Position method of root-finding.<br>
    !>
    !>  \details
    !>  See the root-finding section in the documentation of [pm_mathRoot](@ref pm_mathRoot) for more details about this root-finding method.<br>
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [false](@ref pm_mathRoot::false) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [brent](@ref pm_mathRoot::brent)<br>
    !>  [false](@ref pm_mathRoot::false)<br>
    !>  [secant](@ref pm_mathRoot::secant)<br>
    !>  [halley](@ref pm_mathRoot::halley)<br>
    !>  [newton](@ref pm_mathRoot::newton)<br>
    !>  [ridders](@ref pm_mathRoot::ridders)<br>
    !>  [toms748](@ref pm_mathRoot::toms748)<br>
    !>  [schroder](@ref pm_mathRoot::schroder)<br>
    !>  [bisection](@ref pm_mathRoot::bisection)<br>
    !>  [brent_type](@ref pm_mathRoot::brent_type)<br>
    !>  [false_type](@ref pm_mathRoot::false_type)<br>
    !>  [secant_type](@ref pm_mathRoot::secant_type)<br>
    !>  [halley_type](@ref pm_mathRoot::halley_type)<br>
    !>  [newton_type](@ref pm_mathRoot::newton_type)<br>
    !>  [ridders_type](@ref pm_mathRoot::ridders_type)<br>
    !>  [toms748_type](@ref pm_mathRoot::toms748_type)<br>
    !>  [schroder_type](@ref pm_mathRoot::schroder_type)<br>
    !>  [bisection_type](@ref pm_mathRoot::bisection_type)<br>
    !>  [iteration_type](@ref pm_mathRoot::iteration_type)<br>
    !>  [bracket_type](@ref pm_mathRoot::bracket_type)<br>
    !>  [hybrid_type](@ref pm_mathRoot::hybrid_type)<br>
    !>  [method_type](@ref pm_mathRoot::method_type)<br>
    !>
    !>  \final{false_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(bracket_type) :: false_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [false_type](@ref pm_mathRoot::false_type) that is exclusively used
    !>  to signify the use of False-Position method of root-finding within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  See the root-finding section in the documentation of [pm_mathRoot](@ref pm_mathRoot) for more details about this root-finding method.<br>
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [brent](@ref pm_mathRoot::brent)<br>
    !>  [false](@ref pm_mathRoot::false)<br>
    !>  [secant](@ref pm_mathRoot::secant)<br>
    !>  [halley](@ref pm_mathRoot::halley)<br>
    !>  [newton](@ref pm_mathRoot::newton)<br>
    !>  [ridders](@ref pm_mathRoot::ridders)<br>
    !>  [toms748](@ref pm_mathRoot::toms748)<br>
    !>  [schroder](@ref pm_mathRoot::schroder)<br>
    !>  [bisection](@ref pm_mathRoot::bisection)<br>
    !>  [brent_type](@ref pm_mathRoot::brent_type)<br>
    !>  [false_type](@ref pm_mathRoot::false_type)<br>
    !>  [secant_type](@ref pm_mathRoot::secant_type)<br>
    !>  [halley_type](@ref pm_mathRoot::halley_type)<br>
    !>  [newton_type](@ref pm_mathRoot::newton_type)<br>
    !>  [ridders_type](@ref pm_mathRoot::ridders_type)<br>
    !>  [toms748_type](@ref pm_mathRoot::toms748_type)<br>
    !>  [schroder_type](@ref pm_mathRoot::schroder_type)<br>
    !>  [bisection_type](@ref pm_mathRoot::bisection_type)<br>
    !>  [iteration_type](@ref pm_mathRoot::iteration_type)<br>
    !>  [bracket_type](@ref pm_mathRoot::bracket_type)<br>
    !>  [hybrid_type](@ref pm_mathRoot::hybrid_type)<br>
    !>  [method_type](@ref pm_mathRoot::method_type)<br>
    !>
    !>  \final{false}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(false_type), parameter :: false = false_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: false
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the use of the Secant method of root-finding.<br>
    !>
    !>  \details
    !>  See the root-finding section in the documentation of [pm_mathRoot](@ref pm_mathRoot) for more details about this root-finding method.<br>
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [secant](@ref pm_mathRoot::secant) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [brent](@ref pm_mathRoot::brent)<br>
    !>  [false](@ref pm_mathRoot::false)<br>
    !>  [secant](@ref pm_mathRoot::secant)<br>
    !>  [halley](@ref pm_mathRoot::halley)<br>
    !>  [newton](@ref pm_mathRoot::newton)<br>
    !>  [ridders](@ref pm_mathRoot::ridders)<br>
    !>  [toms748](@ref pm_mathRoot::toms748)<br>
    !>  [schroder](@ref pm_mathRoot::schroder)<br>
    !>  [bisection](@ref pm_mathRoot::bisection)<br>
    !>  [brent_type](@ref pm_mathRoot::brent_type)<br>
    !>  [false_type](@ref pm_mathRoot::false_type)<br>
    !>  [secant_type](@ref pm_mathRoot::secant_type)<br>
    !>  [halley_type](@ref pm_mathRoot::halley_type)<br>
    !>  [newton_type](@ref pm_mathRoot::newton_type)<br>
    !>  [ridders_type](@ref pm_mathRoot::ridders_type)<br>
    !>  [toms748_type](@ref pm_mathRoot::toms748_type)<br>
    !>  [schroder_type](@ref pm_mathRoot::schroder_type)<br>
    !>  [bisection_type](@ref pm_mathRoot::bisection_type)<br>
    !>  [iteration_type](@ref pm_mathRoot::iteration_type)<br>
    !>  [bracket_type](@ref pm_mathRoot::bracket_type)<br>
    !>  [hybrid_type](@ref pm_mathRoot::hybrid_type)<br>
    !>  [method_type](@ref pm_mathRoot::method_type)<br>
    !>
    !>  \final{secant_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(iteration_type) :: secant_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [secant_type](@ref pm_mathRoot::secant_type) that is exclusively used
    !>  to signify the use of Secant method of root-finding within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  See the root-finding section in the documentation of [pm_mathRoot](@ref pm_mathRoot) for more details about this root-finding method.<br>
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [brent](@ref pm_mathRoot::brent)<br>
    !>  [false](@ref pm_mathRoot::false)<br>
    !>  [secant](@ref pm_mathRoot::secant)<br>
    !>  [halley](@ref pm_mathRoot::halley)<br>
    !>  [newton](@ref pm_mathRoot::newton)<br>
    !>  [ridders](@ref pm_mathRoot::ridders)<br>
    !>  [toms748](@ref pm_mathRoot::toms748)<br>
    !>  [schroder](@ref pm_mathRoot::schroder)<br>
    !>  [bisection](@ref pm_mathRoot::bisection)<br>
    !>  [brent_type](@ref pm_mathRoot::brent_type)<br>
    !>  [false_type](@ref pm_mathRoot::false_type)<br>
    !>  [secant_type](@ref pm_mathRoot::secant_type)<br>
    !>  [halley_type](@ref pm_mathRoot::halley_type)<br>
    !>  [newton_type](@ref pm_mathRoot::newton_type)<br>
    !>  [ridders_type](@ref pm_mathRoot::ridders_type)<br>
    !>  [toms748_type](@ref pm_mathRoot::toms748_type)<br>
    !>  [schroder_type](@ref pm_mathRoot::schroder_type)<br>
    !>  [bisection_type](@ref pm_mathRoot::bisection_type)<br>
    !>  [iteration_type](@ref pm_mathRoot::iteration_type)<br>
    !>  [bracket_type](@ref pm_mathRoot::bracket_type)<br>
    !>  [hybrid_type](@ref pm_mathRoot::hybrid_type)<br>
    !>  [method_type](@ref pm_mathRoot::method_type)<br>
    !>
    !>  \final{secant}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(secant_type), parameter :: secant = secant_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: secant
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the use of the Newton method of root-finding.<br>
    !>
    !>  \details
    !>  See the root-finding section in the documentation of [pm_mathRoot](@ref pm_mathRoot) for more details about this root-finding method.<br>
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [newton](@ref pm_mathRoot::newton) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [brent](@ref pm_mathRoot::brent)<br>
    !>  [false](@ref pm_mathRoot::false)<br>
    !>  [secant](@ref pm_mathRoot::secant)<br>
    !>  [halley](@ref pm_mathRoot::halley)<br>
    !>  [newton](@ref pm_mathRoot::newton)<br>
    !>  [ridders](@ref pm_mathRoot::ridders)<br>
    !>  [toms748](@ref pm_mathRoot::toms748)<br>
    !>  [schroder](@ref pm_mathRoot::schroder)<br>
    !>  [bisection](@ref pm_mathRoot::bisection)<br>
    !>  [brent_type](@ref pm_mathRoot::brent_type)<br>
    !>  [false_type](@ref pm_mathRoot::false_type)<br>
    !>  [secant_type](@ref pm_mathRoot::secant_type)<br>
    !>  [halley_type](@ref pm_mathRoot::halley_type)<br>
    !>  [newton_type](@ref pm_mathRoot::newton_type)<br>
    !>  [ridders_type](@ref pm_mathRoot::ridders_type)<br>
    !>  [toms748_type](@ref pm_mathRoot::toms748_type)<br>
    !>  [schroder_type](@ref pm_mathRoot::schroder_type)<br>
    !>  [bisection_type](@ref pm_mathRoot::bisection_type)<br>
    !>  [iteration_type](@ref pm_mathRoot::iteration_type)<br>
    !>  [bracket_type](@ref pm_mathRoot::bracket_type)<br>
    !>  [hybrid_type](@ref pm_mathRoot::hybrid_type)<br>
    !>  [method_type](@ref pm_mathRoot::method_type)<br>
    !>
    !>  \final{newton_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(iteration_type) :: newton_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [newton_type](@ref pm_mathRoot::newton_type) that is exclusively used
    !>  to signify the use of Newton method of root-finding within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  See the root-finding section in the documentation of [pm_mathRoot](@ref pm_mathRoot) for more details about this root-finding method.<br>
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [brent](@ref pm_mathRoot::brent)<br>
    !>  [false](@ref pm_mathRoot::false)<br>
    !>  [secant](@ref pm_mathRoot::secant)<br>
    !>  [halley](@ref pm_mathRoot::halley)<br>
    !>  [newton](@ref pm_mathRoot::newton)<br>
    !>  [ridders](@ref pm_mathRoot::ridders)<br>
    !>  [toms748](@ref pm_mathRoot::toms748)<br>
    !>  [schroder](@ref pm_mathRoot::schroder)<br>
    !>  [bisection](@ref pm_mathRoot::bisection)<br>
    !>  [brent_type](@ref pm_mathRoot::brent_type)<br>
    !>  [false_type](@ref pm_mathRoot::false_type)<br>
    !>  [secant_type](@ref pm_mathRoot::secant_type)<br>
    !>  [halley_type](@ref pm_mathRoot::halley_type)<br>
    !>  [newton_type](@ref pm_mathRoot::newton_type)<br>
    !>  [ridders_type](@ref pm_mathRoot::ridders_type)<br>
    !>  [toms748_type](@ref pm_mathRoot::toms748_type)<br>
    !>  [schroder_type](@ref pm_mathRoot::schroder_type)<br>
    !>  [bisection_type](@ref pm_mathRoot::bisection_type)<br>
    !>  [iteration_type](@ref pm_mathRoot::iteration_type)<br>
    !>  [bracket_type](@ref pm_mathRoot::bracket_type)<br>
    !>  [hybrid_type](@ref pm_mathRoot::hybrid_type)<br>
    !>  [method_type](@ref pm_mathRoot::method_type)<br>
    !>
    !>  \final{newton}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(newton_type), parameter :: newton = newton_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: newton
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the use of the Halley method of root-finding.<br>
    !>
    !>  \details
    !>  See the root-finding section in the documentation of [pm_mathRoot](@ref pm_mathRoot) for more details about this root-finding method.<br>
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [halley](@ref pm_mathRoot::halley) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [brent](@ref pm_mathRoot::brent)<br>
    !>  [false](@ref pm_mathRoot::false)<br>
    !>  [secant](@ref pm_mathRoot::secant)<br>
    !>  [halley](@ref pm_mathRoot::halley)<br>
    !>  [newton](@ref pm_mathRoot::newton)<br>
    !>  [ridders](@ref pm_mathRoot::ridders)<br>
    !>  [toms748](@ref pm_mathRoot::toms748)<br>
    !>  [schroder](@ref pm_mathRoot::schroder)<br>
    !>  [bisection](@ref pm_mathRoot::bisection)<br>
    !>  [brent_type](@ref pm_mathRoot::brent_type)<br>
    !>  [false_type](@ref pm_mathRoot::false_type)<br>
    !>  [secant_type](@ref pm_mathRoot::secant_type)<br>
    !>  [halley_type](@ref pm_mathRoot::halley_type)<br>
    !>  [newton_type](@ref pm_mathRoot::newton_type)<br>
    !>  [ridders_type](@ref pm_mathRoot::ridders_type)<br>
    !>  [toms748_type](@ref pm_mathRoot::toms748_type)<br>
    !>  [schroder_type](@ref pm_mathRoot::schroder_type)<br>
    !>  [bisection_type](@ref pm_mathRoot::bisection_type)<br>
    !>  [iteration_type](@ref pm_mathRoot::iteration_type)<br>
    !>  [bracket_type](@ref pm_mathRoot::bracket_type)<br>
    !>  [hybrid_type](@ref pm_mathRoot::hybrid_type)<br>
    !>  [method_type](@ref pm_mathRoot::method_type)<br>
    !>
    !>  \final{halley_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(iteration_type) :: halley_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [halley_type](@ref pm_mathRoot::halley_type) that is exclusively used
    !>  to signify the use of Halley method of root-finding within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  See the root-finding section in the documentation of [pm_mathRoot](@ref pm_mathRoot) for more details about this root-finding method.<br>
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [brent](@ref pm_mathRoot::brent)<br>
    !>  [false](@ref pm_mathRoot::false)<br>
    !>  [secant](@ref pm_mathRoot::secant)<br>
    !>  [halley](@ref pm_mathRoot::halley)<br>
    !>  [newton](@ref pm_mathRoot::newton)<br>
    !>  [ridders](@ref pm_mathRoot::ridders)<br>
    !>  [toms748](@ref pm_mathRoot::toms748)<br>
    !>  [schroder](@ref pm_mathRoot::schroder)<br>
    !>  [bisection](@ref pm_mathRoot::bisection)<br>
    !>  [brent_type](@ref pm_mathRoot::brent_type)<br>
    !>  [false_type](@ref pm_mathRoot::false_type)<br>
    !>  [secant_type](@ref pm_mathRoot::secant_type)<br>
    !>  [halley_type](@ref pm_mathRoot::halley_type)<br>
    !>  [newton_type](@ref pm_mathRoot::newton_type)<br>
    !>  [ridders_type](@ref pm_mathRoot::ridders_type)<br>
    !>  [toms748_type](@ref pm_mathRoot::toms748_type)<br>
    !>  [schroder_type](@ref pm_mathRoot::schroder_type)<br>
    !>  [bisection_type](@ref pm_mathRoot::bisection_type)<br>
    !>  [iteration_type](@ref pm_mathRoot::iteration_type)<br>
    !>  [bracket_type](@ref pm_mathRoot::bracket_type)<br>
    !>  [hybrid_type](@ref pm_mathRoot::hybrid_type)<br>
    !>  [method_type](@ref pm_mathRoot::method_type)<br>
    !>
    !>  \final{halley}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(halley_type), parameter :: halley = halley_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: halley
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the use of the Schroder method of root-finding.<br>
    !>
    !>  \details
    !>  See the root-finding section in the documentation of [pm_mathRoot](@ref pm_mathRoot) for more details about this root-finding method.<br>
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [schroder](@ref pm_mathRoot::schroder) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [brent](@ref pm_mathRoot::brent)<br>
    !>  [false](@ref pm_mathRoot::false)<br>
    !>  [secant](@ref pm_mathRoot::secant)<br>
    !>  [halley](@ref pm_mathRoot::halley)<br>
    !>  [newton](@ref pm_mathRoot::newton)<br>
    !>  [ridders](@ref pm_mathRoot::ridders)<br>
    !>  [toms748](@ref pm_mathRoot::toms748)<br>
    !>  [schroder](@ref pm_mathRoot::schroder)<br>
    !>  [bisection](@ref pm_mathRoot::bisection)<br>
    !>  [brent_type](@ref pm_mathRoot::brent_type)<br>
    !>  [false_type](@ref pm_mathRoot::false_type)<br>
    !>  [secant_type](@ref pm_mathRoot::secant_type)<br>
    !>  [halley_type](@ref pm_mathRoot::halley_type)<br>
    !>  [newton_type](@ref pm_mathRoot::newton_type)<br>
    !>  [ridders_type](@ref pm_mathRoot::ridders_type)<br>
    !>  [toms748_type](@ref pm_mathRoot::toms748_type)<br>
    !>  [schroder_type](@ref pm_mathRoot::schroder_type)<br>
    !>  [bisection_type](@ref pm_mathRoot::bisection_type)<br>
    !>  [iteration_type](@ref pm_mathRoot::iteration_type)<br>
    !>  [bracket_type](@ref pm_mathRoot::bracket_type)<br>
    !>  [hybrid_type](@ref pm_mathRoot::hybrid_type)<br>
    !>  [method_type](@ref pm_mathRoot::method_type)<br>
    !>
    !>  \final{schroder_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(iteration_type) :: schroder_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [schroder_type](@ref pm_mathRoot::schroder_type) that is exclusively used
    !>  to signify the use of Schroder method of root-finding within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  See the root-finding section in the documentation of [pm_mathRoot](@ref pm_mathRoot) for more details about this root-finding method.<br>
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [brent](@ref pm_mathRoot::brent)<br>
    !>  [false](@ref pm_mathRoot::false)<br>
    !>  [secant](@ref pm_mathRoot::secant)<br>
    !>  [halley](@ref pm_mathRoot::halley)<br>
    !>  [newton](@ref pm_mathRoot::newton)<br>
    !>  [ridders](@ref pm_mathRoot::ridders)<br>
    !>  [toms748](@ref pm_mathRoot::toms748)<br>
    !>  [schroder](@ref pm_mathRoot::schroder)<br>
    !>  [bisection](@ref pm_mathRoot::bisection)<br>
    !>  [brent_type](@ref pm_mathRoot::brent_type)<br>
    !>  [false_type](@ref pm_mathRoot::false_type)<br>
    !>  [secant_type](@ref pm_mathRoot::secant_type)<br>
    !>  [halley_type](@ref pm_mathRoot::halley_type)<br>
    !>  [newton_type](@ref pm_mathRoot::newton_type)<br>
    !>  [ridders_type](@ref pm_mathRoot::ridders_type)<br>
    !>  [toms748_type](@ref pm_mathRoot::toms748_type)<br>
    !>  [schroder_type](@ref pm_mathRoot::schroder_type)<br>
    !>  [bisection_type](@ref pm_mathRoot::bisection_type)<br>
    !>  [iteration_type](@ref pm_mathRoot::iteration_type)<br>
    !>  [bracket_type](@ref pm_mathRoot::bracket_type)<br>
    !>  [hybrid_type](@ref pm_mathRoot::hybrid_type)<br>
    !>  [method_type](@ref pm_mathRoot::method_type)<br>
    !>
    !>  \final{schroder}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(schroder_type), parameter :: schroder = schroder_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: schroder
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the use of the Ridders method of root-finding.<br>
    !>
    !>  \details
    !>  See the root-finding section in the documentation of [pm_mathRoot](@ref pm_mathRoot) for more details about this root-finding method.<br>
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [ridders](@ref pm_mathRoot::ridders) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [brent](@ref pm_mathRoot::brent)<br>
    !>  [false](@ref pm_mathRoot::false)<br>
    !>  [secant](@ref pm_mathRoot::secant)<br>
    !>  [halley](@ref pm_mathRoot::halley)<br>
    !>  [newton](@ref pm_mathRoot::newton)<br>
    !>  [ridders](@ref pm_mathRoot::ridders)<br>
    !>  [toms748](@ref pm_mathRoot::toms748)<br>
    !>  [schroder](@ref pm_mathRoot::schroder)<br>
    !>  [bisection](@ref pm_mathRoot::bisection)<br>
    !>  [brent_type](@ref pm_mathRoot::brent_type)<br>
    !>  [false_type](@ref pm_mathRoot::false_type)<br>
    !>  [secant_type](@ref pm_mathRoot::secant_type)<br>
    !>  [halley_type](@ref pm_mathRoot::halley_type)<br>
    !>  [newton_type](@ref pm_mathRoot::newton_type)<br>
    !>  [ridders_type](@ref pm_mathRoot::ridders_type)<br>
    !>  [toms748_type](@ref pm_mathRoot::toms748_type)<br>
    !>  [schroder_type](@ref pm_mathRoot::schroder_type)<br>
    !>  [bisection_type](@ref pm_mathRoot::bisection_type)<br>
    !>  [iteration_type](@ref pm_mathRoot::iteration_type)<br>
    !>  [bracket_type](@ref pm_mathRoot::bracket_type)<br>
    !>  [hybrid_type](@ref pm_mathRoot::hybrid_type)<br>
    !>  [method_type](@ref pm_mathRoot::method_type)<br>
    !>
    !>  \final{ridders_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(hybrid_type) :: ridders_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [ridders_type](@ref pm_mathRoot::ridders_type) that is exclusively used
    !>  to signify the use of Ridders method of root-finding within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  See the root-finding section in the documentation of [pm_mathRoot](@ref pm_mathRoot) for more details about this root-finding method.<br>
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [brent](@ref pm_mathRoot::brent)<br>
    !>  [false](@ref pm_mathRoot::false)<br>
    !>  [secant](@ref pm_mathRoot::secant)<br>
    !>  [halley](@ref pm_mathRoot::halley)<br>
    !>  [newton](@ref pm_mathRoot::newton)<br>
    !>  [ridders](@ref pm_mathRoot::ridders)<br>
    !>  [toms748](@ref pm_mathRoot::toms748)<br>
    !>  [schroder](@ref pm_mathRoot::schroder)<br>
    !>  [bisection](@ref pm_mathRoot::bisection)<br>
    !>  [brent_type](@ref pm_mathRoot::brent_type)<br>
    !>  [false_type](@ref pm_mathRoot::false_type)<br>
    !>  [secant_type](@ref pm_mathRoot::secant_type)<br>
    !>  [halley_type](@ref pm_mathRoot::halley_type)<br>
    !>  [newton_type](@ref pm_mathRoot::newton_type)<br>
    !>  [ridders_type](@ref pm_mathRoot::ridders_type)<br>
    !>  [toms748_type](@ref pm_mathRoot::toms748_type)<br>
    !>  [schroder_type](@ref pm_mathRoot::schroder_type)<br>
    !>  [bisection_type](@ref pm_mathRoot::bisection_type)<br>
    !>  [iteration_type](@ref pm_mathRoot::iteration_type)<br>
    !>  [bracket_type](@ref pm_mathRoot::bracket_type)<br>
    !>  [hybrid_type](@ref pm_mathRoot::hybrid_type)<br>
    !>  [method_type](@ref pm_mathRoot::method_type)<br>
    !>
    !>  \final{ridders}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(ridders_type), parameter :: ridders = ridders_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: ridders
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the use of the Bisection method of root-finding.<br>
    !>
    !>  \details
    !>  See the root-finding section in the documentation of [pm_mathRoot](@ref pm_mathRoot) for more details about this root-finding method.<br>
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [bisection](@ref pm_mathRoot::bisection) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [brent](@ref pm_mathRoot::brent)<br>
    !>  [false](@ref pm_mathRoot::false)<br>
    !>  [secant](@ref pm_mathRoot::secant)<br>
    !>  [halley](@ref pm_mathRoot::halley)<br>
    !>  [newton](@ref pm_mathRoot::newton)<br>
    !>  [ridders](@ref pm_mathRoot::ridders)<br>
    !>  [toms748](@ref pm_mathRoot::toms748)<br>
    !>  [schroder](@ref pm_mathRoot::schroder)<br>
    !>  [bisection](@ref pm_mathRoot::bisection)<br>
    !>  [brent_type](@ref pm_mathRoot::brent_type)<br>
    !>  [false_type](@ref pm_mathRoot::false_type)<br>
    !>  [secant_type](@ref pm_mathRoot::secant_type)<br>
    !>  [halley_type](@ref pm_mathRoot::halley_type)<br>
    !>  [newton_type](@ref pm_mathRoot::newton_type)<br>
    !>  [ridders_type](@ref pm_mathRoot::ridders_type)<br>
    !>  [toms748_type](@ref pm_mathRoot::toms748_type)<br>
    !>  [schroder_type](@ref pm_mathRoot::schroder_type)<br>
    !>  [bisection_type](@ref pm_mathRoot::bisection_type)<br>
    !>  [iteration_type](@ref pm_mathRoot::iteration_type)<br>
    !>  [bracket_type](@ref pm_mathRoot::bracket_type)<br>
    !>  [hybrid_type](@ref pm_mathRoot::hybrid_type)<br>
    !>  [method_type](@ref pm_mathRoot::method_type)<br>
    !>
    !>  \final{bisection_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(bracket_type) :: bisection_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [bisection_type](@ref pm_mathRoot::bisection_type) that is exclusively used
    !>  to signify the use of Bisection method of root-finding within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  See the root-finding section in the documentation of [pm_mathRoot](@ref pm_mathRoot) for more details about this root-finding method.<br>
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [brent](@ref pm_mathRoot::brent)<br>
    !>  [false](@ref pm_mathRoot::false)<br>
    !>  [secant](@ref pm_mathRoot::secant)<br>
    !>  [halley](@ref pm_mathRoot::halley)<br>
    !>  [newton](@ref pm_mathRoot::newton)<br>
    !>  [ridders](@ref pm_mathRoot::ridders)<br>
    !>  [toms748](@ref pm_mathRoot::toms748)<br>
    !>  [schroder](@ref pm_mathRoot::schroder)<br>
    !>  [bisection](@ref pm_mathRoot::bisection)<br>
    !>  [brent_type](@ref pm_mathRoot::brent_type)<br>
    !>  [false_type](@ref pm_mathRoot::false_type)<br>
    !>  [secant_type](@ref pm_mathRoot::secant_type)<br>
    !>  [halley_type](@ref pm_mathRoot::halley_type)<br>
    !>  [newton_type](@ref pm_mathRoot::newton_type)<br>
    !>  [ridders_type](@ref pm_mathRoot::ridders_type)<br>
    !>  [toms748_type](@ref pm_mathRoot::toms748_type)<br>
    !>  [schroder_type](@ref pm_mathRoot::schroder_type)<br>
    !>  [bisection_type](@ref pm_mathRoot::bisection_type)<br>
    !>  [iteration_type](@ref pm_mathRoot::iteration_type)<br>
    !>  [bracket_type](@ref pm_mathRoot::bracket_type)<br>
    !>  [hybrid_type](@ref pm_mathRoot::hybrid_type)<br>
    !>  [method_type](@ref pm_mathRoot::method_type)<br>
    !>
    !>  \final{bisection}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(bisection_type), parameter :: bisection = bisection_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: bisection
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a root of a specified continuous real-valued one-dimensional
    !>  mathematical function such that \f$f(\mathrm{root}) = 0\f$ with the user-specified or the default root-finding method.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_mathRoot](@ref pm_mathRoot) for details of the various root-finding methods.<br>
    !>
    !>  \note
    !>  <b>Which root-finding method should I use among all available?</b><br>
    !>  The [brent](@ref pm_mathRoot::brent) method is the recommended **algorithm of choice**
    !>  for general one-dimensional root-finding problems **without** knowledge of function **derivative** (gradient).<br>
    !>  The [Newton-Raphson method](@ref pm_mathRoot::newton) is the recommended **algorithm of choice**
    !>  for general one-dimensional root-finding problems **with** knowledge of function **derivative** (gradient).<br>
    !>
    !>  \param[in]      method  :   The input scalar constant that can be one of the following:<br>
    !>                              <ol>
    !>                                  <li>    The constant [false](@ref pm_mathRoot::false) or an object of type [false_type](@ref pm_mathRoot::false_type),
    !>                                          signifying the use of the **False-Position** root-finding method within the algorithm.<br>
    !>                                  <li>    The constant [bisection](@ref pm_mathRoot::bisection) or an object of type [bisection_type](@ref pm_mathRoot::bisection_type),
    !>                                          signifying the use of the **Bisection** root-finding method within the algorithm.<br>
    !>                                  <li>    The constant [secant](@ref pm_mathRoot::secant) or an object of type [secant_type](@ref pm_mathRoot::secant_type),
    !>                                          signifying the use of the **Secant** root-finding method within the algorithm.<br>
    !>                                  <li>    The constant [brent](@ref pm_mathRoot::brent) or an object of type [brent_type](@ref pm_mathRoot::brent_type),
    !>                                          signifying the use of the **Brent** root-finding method within the algorithm.<br>
    !>                                  <li>    The constant [ridders](@ref pm_mathRoot::ridders) or an object of type [ridders_type](@ref pm_mathRoot::ridders_type),
    !>                                          signifying the use of the **Ridders** root-finding method within the algorithm.<br>
    !>                                  <li>    The constant [toms748](@ref pm_mathRoot::toms748) or an object of type [toms748_type](@ref pm_mathRoot::toms748_type),
    !>                                          signifying the use of the **TOMS748** root-finding method within the algorithm.<br>
    !>                                  <li>    The constant [newton](@ref pm_mathRoot::newton) or an object of type [newton_type](@ref pm_mathRoot::newton_type),
    !>                                          signifying the use of the **Newton** root-finding method within the algorithm.<br>
    !>                                  <li>    The constant [halley](@ref pm_mathRoot::halley) or an object of type [halley_type](@ref pm_mathRoot::halley_type),
    !>                                          signifying the use of the **Halley** root-finding method within the algorithm.<br>
    !>                                  <li>    The constant [schroder](@ref pm_mathRoot::schroder) or an object of type [schroder_type](@ref pm_mathRoot::schroder_type),
    !>                                          signifying the use of the **Schroder** root-finding method within the algorithm.<br>
    !>                              </ol>
    !>                              (**optional**, default = [brent](@ref pm_mathRoot::brent). Note the implicit constraint this default option sets on the input `getFunc()` interface below.)
    !>  \param[in]      getFunc :   The `external` user-specified function whose interface depends on the specified value for the input argument `method`.<br>
    !>                              <ol>
    !>                                  <li>    If the specified `method` is any of the following,<br>
    !>                                          <ol>
    !>                                              <li>    the constant [brent](@ref pm_mathRoot::brent) or an object of type [brent_type](@ref pm_mathRoot::brent_type),<br>
    !>                                              <li>    the constant [false](@ref pm_mathRoot::false) or an object of type [false_type](@ref pm_mathRoot::false_type),<br>
    !>                                              <li>    the constant [secant](@ref pm_mathRoot::secant) or an object of type [secant_type](@ref pm_mathRoot::secant_type),<br>
    !>                                              <li>    the constant [ridders](@ref pm_mathRoot::ridders) or an object of type [ridders_type](@ref pm_mathRoot::ridders_type),<br>
    !>                                              <li>    the constant [toms748](@ref pm_mathRoot::toms748) or an object of type [toms748_type](@ref pm_mathRoot::toms748_type),<br>
    !>                                              <li>    the constant [bisection](@ref pm_mathRoot::bisection) or an object of type [bisection_type](@ref pm_mathRoot::bisection_type),<br>
    !>                                          </ol>
    !>                                          none of which require the derivative of the target function to find the function root,
    !>                                          then `getFunc()` must take a single scalar input argument `x` of the same type and kind as the output argument `root` (below).<br>
    !>                                          On output, `getFunc()` must return a scalar of the same type and kind as the function input argument `x`,
    !>                                          containing the value of the target function evaluated at the specified input `x`.<br>
    !>                                          The following illustrates the generic interface of `getFunc()` for the above values of `method`,
    !>                                          \code{.F90}
    !>                                              function getFunc(x) result(func)
    !>                                                  real(RKG)   , intent(in)    :: x
    !>                                                  real(RKG)                   :: func
    !>                                              end function
    !>                                          \endcode
    !>                                          where `RKG` refers to the kind of the output argument `root`.<br>
    !>                                  <li>    If the specified `method` is any of the following,
    !>                                          <ol>
    !>                                              <li>    the constant [newton](@ref pm_mathRoot::newton) or an object of type [newton_type](@ref pm_mathRoot::newton_type),<br>
    !>                                              <li>    the constant [halley](@ref pm_mathRoot::halley) or an object of type [halley_type](@ref pm_mathRoot::halley_type),<br>
    !>                                              <li>    the constant [schroder](@ref pm_mathRoot::schroder) or an object of type [schroder_type](@ref pm_mathRoot::schroder_type),<br>
    !>                                          </ol>
    !>                                          all of which either require the first (Newton) or higher derivatives (Halley/Schroder) of the target function to find the function root,
    !>                                          then `getFunc()` must take two arguments:<br>
    !>                                          <ol>
    !>                                              <li>    A scalar input argument `x` of the same type and kind as the output argument `root` (below).<br>
    !>                                              <li>    A scalar input argument `order` of type `integer` of default kind \IK.<br>
    !>                                          </ol>
    !>                                          On output, `getFunc()` must return a scalar of the same type and kind as the function input argument `x`,
    !>                                          containing the derivative of the specified input `order` of the target function evaluated at the specified input `x`:<br>
    !>                                          <ol>
    !>                                              <li>    An input `order` value of `0` must yield the target function value at the input `x`.<br>
    !>                                              <li>    An input `order` value of `1` must yield the target function first derivative at the input `x`.<br>
    !>                                              <li>    An input `order` value of `2` must yield the target function second derivative at the input `x`.<br>
    !>                                          </ol>
    !>                                          The following illustrates the generic interface of `getFunc()` for the above values of `method`,
    !>                                          \code{.F90}
    !>                                              function getFunc(x, order) result(func)
    !>                                                  real(RKG)   , intent(in)    :: x
    !>                                                  integer(IK) , intent(in)    :: order
    !>                                                  real(RKG)                   :: func
    !>                                              end function
    !>                                          \endcode
    !>                                          where `RKG` refers to the kind of the output argument `root`.<br>
    !>                              </ol>
    !>  \param[inout]   root    :   The input/output scalar of type `real` of kind \RKALL.<br>
    !>                              On input,<br>
    !>                              <ol>
    !>                                  <li>    If the specified `method` is of type [newton_type](@ref pm_mathRoot::newton_type) or [halley_type](@ref pm_mathRoot::halley_type),
    !>                                          then `root` must contain the best guess initial value for the function `root` within the bracket specified by the input arguments `lb` and `ub`.<br>
    !>                                  <li>    Otherwise, any input value for `root` is ignored for all other root-finding methods.<br>
    !>                              </ol>
    !>                              On output, `root` contains the best approximation to the root of the user-specified
    !>                              function `getFunc()` within the specified search interval `[lb, ub]` and tolerance threshold `abstol`.<br>
    !>  \param[in]      lb      :   The input scalar of same type and kind as the output `root`, representing the lower bound of the bracket (search) interval.<br>
    !>  \param[in]      ub      :   The input scalar of same type and kind as the output `root`, representing the upper bound of the bracket (search) interval.<br>
    !>  \param[in]      abstol  :   The input scalar of same type and kind as the output `root`, representing the absolute tolerance used as the stopping criterion of the search.<br>
    !>                              The iterations of the specified input `method` continue until the search interval becomes smaller than `abstol` in absolute units.<br>
    !>                              Care must be taken for specifying a reasonable value for `abstol` (see the warnings below).<br>
    !>                              If no suitable value for `abstol` is known a priori, try `abstol = epsilon(0._RKG)**.8 * (abs(lb) + abs(ub))`
    !>                              where `RKG` refers to the kind of the output argument `root`.<br>
    !>                              (**optional**, default = `epsilon(0._RKG)**.8 * (abs(lb) + abs(ub))`)
    !>  \param[out]     neval   :   The output scalar argument of type `integer` of default kind \IK, containing the total number of `getFunc()` function calls made by the algorithm.<br>
    !>                              <ol>
    !>                                  <li>    A positive output `neval` implies the **successful** convergence of the algorithm to the function root after `+neval` function evaluations.<br>
    !>                                  <li>    A negative output `neval` implies the **failure** of convergence of the algorithm to the function root after `-neval` function evaluations.<br>
    !>                                  <li>    A zero output `neval` only occurs if the root lies at the search bracket boundaries specified by the input argument `lb` or `ub`.<br>
    !>                              </ol>
    !>                              Convergence failures rarely occur. If they do, setting the input arguments `abstol` and/or `niter` to larger values may resolve the failure.<br>
    !>  \param[in]      niter   :   The input scalar of type `integer` of default kind \IK, representing the **maximum number of iterations allowed** for the specified `method` to find the root.<br>
    !>                              The default number of steps `niter` within the algorithm is a compile-time constant that depends on the `kind` of the `real` input arguments.<br>
    !>                              (**optional**, default = `ceiling(2 * precision(lb) / log10(2.))`)
    !>
    !>  \interface{getRoot}
    !>  \code{.F90}
    !>
    !>      use pm_mathRoot, only: getRoot
    !>
    !>      root = getRoot(getFunc, lb, ub, abstol = abstol, neval = neval, niter = niter)
    !>      root = getRoot(method, getFunc, lb, ub, abstol = abstol, neval = neval, niter = niter)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `lb < ub` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < abstol` must hold for the corresponding input arguments.<br>
    !>  The condition `abstol < (ub - lb)` must hold for the corresponding input arguments.<br>
    !>  The condition `abs(getFunc(lb), lf) < abstol` must hold for the corresponding input arguments.<br>
    !>  The condition `abs(getFunc(ub), uf) < abstol` must hold for the corresponding input arguments.<br>
    !>  The condition `sign(lf, 1.) /= sign(uf, 1.)` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < niter` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warning
    !>  It is crucial to keep in mind that computers use a fixed number of binary digits to represent floating-point numbers.<br>
    !>  While the user-specified function might analytically pass through zero, it is possible that its computed value is never zero, for any floating-point argument.<br>
    !>  One must also decide on what accuracy on the root is attainable.<br>
    !>  For example, convergence to within \f$10^{−6}\f$ in absolute value is reasonable when the root lies near \f$1\f$,
    !>  but unachievable if the root lies near an extremely large number near floating-point overflow.<br>
    !>
    !>  \impure
    !>
    !>  \recursive
    !>
    !>  \see
    !>  [getRoot](@ref pm_mathRoot::getRoot)<br>
    !>  [setRoot](@ref pm_mathRoot::setRoot)<br>
    !>  [getPolyRoot](@ref pm_polynomial::getPolyRoot)<br>
    !>  [setPolyRoot](@ref pm_polynomial::setPolyRoot)<br>
    !>
    !>  \example{getRoot}
    !>  \include{lineno} example/pm_mathRoot/getRoot/main.F90
    !>  \compilef{getRoot}
    !>  \output{getRoot}
    !>  \include{lineno} example/pm_mathRoot/getRoot/main.out.F90
    !>
    !>  \test
    !>  [test_pm_mathRoot](@ref test_pm_mathRoot)
    !>
    !>  \final{getRoot}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

    ! Default

    interface getRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module function getRootDef_RK5(getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootDef_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK4_ENABLED
    recursive module function getRootDef_RK4(getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootDef_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK3_ENABLED
    recursive module function getRootDef_RK3(getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootDef_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK2_ENABLED
    recursive module function getRootDef_RK2(getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootDef_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK1_ENABLED
    recursive module function getRootDef_RK1(getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootDef_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! False

    interface getRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module function getRootFalse_RK5(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootFalse_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(false_type)        , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK4_ENABLED
    recursive module function getRootFalse_RK4(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootFalse_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(false_type)        , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK3_ENABLED
    recursive module function getRootFalse_RK3(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootFalse_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(false_type)        , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK2_ENABLED
    recursive module function getRootFalse_RK2(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootFalse_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(false_type)        , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK1_ENABLED
    recursive module function getRootFalse_RK1(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootFalse_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(false_type)        , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Bisection

    interface getRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module function getRootBisection_RK5(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootBisection_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(bisection_type)    , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK4_ENABLED
    recursive module function getRootBisection_RK4(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootBisection_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(bisection_type)    , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK3_ENABLED
    recursive module function getRootBisection_RK3(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootBisection_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(bisection_type)    , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK2_ENABLED
    recursive module function getRootBisection_RK2(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootBisection_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(bisection_type)    , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK1_ENABLED
    recursive module function getRootBisection_RK1(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootBisection_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(bisection_type)    , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Secant

    interface getRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module function getRootSecant_RK5(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootSecant_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(secant_type)       , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK4_ENABLED
    recursive module function getRootSecant_RK4(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootSecant_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(secant_type)       , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK3_ENABLED
    recursive module function getRootSecant_RK3(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootSecant_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(secant_type)       , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK2_ENABLED
    recursive module function getRootSecant_RK2(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootSecant_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(secant_type)       , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK1_ENABLED
    recursive module function getRootSecant_RK1(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootSecant_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(secant_type)       , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Brent

    interface getRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module function getRootBrent_RK5(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootBrent_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(brent_type)        , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK4_ENABLED
    recursive module function getRootBrent_RK4(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootBrent_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(brent_type)        , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK3_ENABLED
    recursive module function getRootBrent_RK3(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootBrent_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(brent_type)        , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK2_ENABLED
    recursive module function getRootBrent_RK2(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootBrent_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(brent_type)        , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK1_ENABLED
    recursive module function getRootBrent_RK1(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootBrent_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(brent_type)        , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Ridders

    interface getRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module function getRootRidders_RK5(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootRidders_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(ridders_type)      , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK4_ENABLED
    recursive module function getRootRidders_RK4(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootRidders_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(ridders_type)      , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK3_ENABLED
    recursive module function getRootRidders_RK3(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootRidders_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(ridders_type)      , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK2_ENABLED
    recursive module function getRootRidders_RK2(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootRidders_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(ridders_type)      , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK1_ENABLED
    recursive module function getRootRidders_RK1(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootRidders_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(ridders_type)      , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! TOMS748

    interface getRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module function getRootTOMS748_RK5(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootTOMS748_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(toms748_type)      , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK4_ENABLED
    recursive module function getRootTOMS748_RK4(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootTOMS748_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(toms748_type)      , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK3_ENABLED
    recursive module function getRootTOMS748_RK3(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootTOMS748_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(toms748_type)      , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK2_ENABLED
    recursive module function getRootTOMS748_RK2(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootTOMS748_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(toms748_type)      , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

#if RK1_ENABLED
    recursive module function getRootTOMS748_RK1(method, getFunc, lb, ub, abstol, neval, niter) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootTOMS748_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(toms748_type)      , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)                                           :: root
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Newton

    interface getRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module function getRootNewton_RK5(method, getFunc, lb, ub, abstol, neval, niter, init) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootNewton_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(newton_type)       , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)               , intent(in)    , optional  :: init
        real(RKG)                                           :: root
    end function
#endif

#if RK4_ENABLED
    recursive module function getRootNewton_RK4(method, getFunc, lb, ub, abstol, neval, niter, init) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootNewton_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(newton_type)       , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)               , intent(in)    , optional  :: init
        real(RKG)                                           :: root
    end function
#endif

#if RK3_ENABLED
    recursive module function getRootNewton_RK3(method, getFunc, lb, ub, abstol, neval, niter, init) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootNewton_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(newton_type)       , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)               , intent(in)    , optional  :: init
        real(RKG)                                           :: root
    end function
#endif

#if RK2_ENABLED
    recursive module function getRootNewton_RK2(method, getFunc, lb, ub, abstol, neval, niter, init) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootNewton_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(newton_type)       , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)               , intent(in)    , optional  :: init
        real(RKG)                                           :: root
    end function
#endif

#if RK1_ENABLED
    recursive module function getRootNewton_RK1(method, getFunc, lb, ub, abstol, neval, niter, init) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootNewton_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(newton_type)       , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)               , intent(in)    , optional  :: init
        real(RKG)                                           :: root
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Halley

    interface getRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module function getRootHalley_RK5(method, getFunc, lb, ub, abstol, neval, niter, init) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootHalley_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(halley_type)       , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)               , intent(in)    , optional  :: init
        real(RKG)                                           :: root
    end function
#endif

#if RK4_ENABLED
    recursive module function getRootHalley_RK4(method, getFunc, lb, ub, abstol, neval, niter, init) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootHalley_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(halley_type)       , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)               , intent(in)    , optional  :: init
        real(RKG)                                           :: root
    end function
#endif

#if RK3_ENABLED
    recursive module function getRootHalley_RK3(method, getFunc, lb, ub, abstol, neval, niter, init) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootHalley_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(halley_type)       , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)               , intent(in)    , optional  :: init
        real(RKG)                                           :: root
    end function
#endif

#if RK2_ENABLED
    recursive module function getRootHalley_RK2(method, getFunc, lb, ub, abstol, neval, niter, init) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootHalley_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(halley_type)       , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)               , intent(in)    , optional  :: init
        real(RKG)                                           :: root
    end function
#endif

#if RK1_ENABLED
    recursive module function getRootHalley_RK1(method, getFunc, lb, ub, abstol, neval, niter, init) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootHalley_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(halley_type)       , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)               , intent(in)    , optional  :: init
        real(RKG)                                           :: root
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Schroder

    interface getRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module function getRootSchroder_RK5(method, getFunc, lb, ub, abstol, neval, niter, init) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootSchroder_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(schroder_type)     , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)               , intent(in)    , optional  :: init
        real(RKG)                                           :: root
    end function
#endif

#if RK4_ENABLED
    recursive module function getRootSchroder_RK4(method, getFunc, lb, ub, abstol, neval, niter, init) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootSchroder_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(schroder_type)     , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)               , intent(in)    , optional  :: init
        real(RKG)                                           :: root
    end function
#endif

#if RK3_ENABLED
    recursive module function getRootSchroder_RK3(method, getFunc, lb, ub, abstol, neval, niter, init) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootSchroder_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(schroder_type)     , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)               , intent(in)    , optional  :: init
        real(RKG)                                           :: root
    end function
#endif

#if RK2_ENABLED
    recursive module function getRootSchroder_RK2(method, getFunc, lb, ub, abstol, neval, niter, init) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootSchroder_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(schroder_type)     , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)               , intent(in)    , optional  :: init
        real(RKG)                                           :: root
    end function
#endif

#if RK1_ENABLED
    recursive module function getRootSchroder_RK1(method, getFunc, lb, ub, abstol, neval, niter, init) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRootSchroder_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(schroder_type)     , intent(in)                :: method
        real(RKG)               , intent(in)                :: lb, ub
        real(RKG)               , intent(in)    , optional  :: abstol
        integer(IK)             , intent(in)    , optional  :: niter
        integer(IK)             , intent(out)   , optional  :: neval
        real(RKG)               , intent(in)    , optional  :: init
        real(RKG)                                           :: root
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a root of a specified continuous real-valued one-dimensional
    !>  mathematical function such that \f$f(\mathrm{root}) = 0\f$ with the user-specified or the default root-finding method.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_mathRoot](@ref pm_mathRoot) for details of the various root-finding methods.<br>
    !>
    !>  \note
    !>  <b>Which root-finding method should I use among all available?</b><br>
    !>  The [brent](@ref pm_mathRoot::brent) method is the recommended **algorithm of choice**
    !>  for general one-dimensional root-finding problems **without** knowledge of function **derivative** (gradient).<br>
    !>  The [Newton-Raphson method](@ref pm_mathRoot::newton) is the recommended **algorithm of choice**
    !>  for general one-dimensional root-finding problems **with** knowledge of function **derivative** (gradient).<br>
    !>
    !>  \param[in]      method  :   The input scalar constant that can be one of the following:<br>
    !>                              <ol>
    !>                                  <li>    The constant [false](@ref pm_mathRoot::false) or an object of type [false_type](@ref pm_mathRoot::false_type),
    !>                                          signifying the use of the **False-Position** root-finding method within the algorithm.<br>
    !>                                  <li>    The constant [bisection](@ref pm_mathRoot::bisection) or an object of type [bisection_type](@ref pm_mathRoot::bisection_type),
    !>                                          signifying the use of the **Bisection** root-finding method within the algorithm.<br>
    !>                                  <li>    The constant [secant](@ref pm_mathRoot::secant) or an object of type [secant_type](@ref pm_mathRoot::secant_type),
    !>                                          signifying the use of the **Secant** root-finding method within the algorithm.<br>
    !>                                  <li>    The constant [brent](@ref pm_mathRoot::brent) or an object of type [brent_type](@ref pm_mathRoot::brent_type),
    !>                                          signifying the use of the **Brent** root-finding method within the algorithm.<br>
    !>                                  <li>    The constant [ridders](@ref pm_mathRoot::ridders) or an object of type [ridders_type](@ref pm_mathRoot::ridders_type),
    !>                                          signifying the use of the **Ridders** root-finding method within the algorithm.<br>
    !>                                  <li>    The constant [toms748](@ref pm_mathRoot::toms748) or an object of type [toms748_type](@ref pm_mathRoot::toms748_type),
    !>                                          signifying the use of the **TOMS748** root-finding method within the algorithm.<br>
    !>                                  <li>    The constant [newton](@ref pm_mathRoot::newton) or an object of type [newton_type](@ref pm_mathRoot::newton_type),
    !>                                          signifying the use of the **Newton** root-finding method within the algorithm.<br>
    !>                                  <li>    The constant [halley](@ref pm_mathRoot::halley) or an object of type [halley_type](@ref pm_mathRoot::halley_type),
    !>                                          signifying the use of the **Halley** root-finding method within the algorithm.<br>
    !>                                  <li>    The constant [schroder](@ref pm_mathRoot::schroder) or an object of type [schroder_type](@ref pm_mathRoot::schroder_type),
    !>                                          signifying the use of the **Schroder** root-finding method within the algorithm.<br>
    !>                              </ol>
    !>  \param[in]      getFunc :   The `external` user-specified function whose interface depends on the specified value for the input argument `method`.<br>
    !>                              <ol>
    !>                                  <li>    If the specified `method` is any of the following,<br>
    !>                                          <ol>
    !>                                              <li>    the constant [brent](@ref pm_mathRoot::brent) or an object of type [brent_type](@ref pm_mathRoot::brent_type),<br>
    !>                                              <li>    the constant [false](@ref pm_mathRoot::false) or an object of type [false_type](@ref pm_mathRoot::false_type),<br>
    !>                                              <li>    the constant [secant](@ref pm_mathRoot::secant) or an object of type [secant_type](@ref pm_mathRoot::secant_type),<br>
    !>                                              <li>    the constant [ridders](@ref pm_mathRoot::ridders) or an object of type [ridders_type](@ref pm_mathRoot::ridders_type),<br>
    !>                                              <li>    the constant [toms748](@ref pm_mathRoot::toms748) or an object of type [toms748_type](@ref pm_mathRoot::toms748_type),<br>
    !>                                              <li>    the constant [bisection](@ref pm_mathRoot::bisection) or an object of type [bisection_type](@ref pm_mathRoot::bisection_type),<br>
    !>                                          </ol>
    !>                                          none of which require the derivative of the target function to find the function root,
    !>                                          then `getFunc()` must take a single scalar input argument `x` of the same type and kind as the output argument `root` (below).<br>
    !>                                          On output, `getFunc()` must return a scalar of the same type and kind as the function input argument `x`,
    !>                                          containing the value of the target function evaluated at the specified input `x`.<br>
    !>                                          The following illustrates the generic interface of `getFunc()` for the above values of `method`,
    !>                                          \code{.F90}
    !>                                              function getFunc(x) result(func)
    !>                                                  real(RKG)   , intent(in)    :: x
    !>                                                  real(RKG)                   :: func
    !>                                              end function
    !>                                          \endcode
    !>                                          where `RKG` refers to the kind of the output argument `root`.<br>
    !>                                  <li>    If the specified `method` is any of the following,
    !>                                          <ol>
    !>                                              <li>    the constant [newton](@ref pm_mathRoot::newton) or an object of type [newton_type](@ref pm_mathRoot::newton_type),<br>
    !>                                              <li>    the constant [halley](@ref pm_mathRoot::halley) or an object of type [halley_type](@ref pm_mathRoot::halley_type),<br>
    !>                                              <li>    the constant [schroder](@ref pm_mathRoot::schroder) or an object of type [schroder_type](@ref pm_mathRoot::schroder_type),<br>
    !>                                          </ol>
    !>                                          all of which either require the first (Newton) or higher derivatives (Halley/Schroder) of the target function to find the function root,
    !>                                          then `getFunc()` must take two arguments:<br>
    !>                                          <ol>
    !>                                              <li>    A scalar input argument `x` of the same type and kind as the output argument `root` (below).<br>
    !>                                              <li>    A scalar input argument `order` of type `integer` of default kind \IK.<br>
    !>                                          </ol>
    !>                                          On output, `getFunc()` must return a scalar of the same type and kind as the function input argument `x`,
    !>                                          containing the derivative of the specified input `order` of the target function evaluated at the specified input `x`:<br>
    !>                                          <ol>
    !>                                              <li>    An input `order` value of `0` must yield the target function value at the input `x`.<br>
    !>                                              <li>    An input `order` value of `1` must yield the target function first derivative at the input `x`.<br>
    !>                                              <li>    An input `order` value of `2` must yield the target function second derivative at the input `x`.<br>
    !>                                          </ol>
    !>                                          The following illustrates the generic interface of `getFunc()` for the above values of `method`,
    !>                                          \code{.F90}
    !>                                              function getFunc(x, order) result(func)
    !>                                                  real(RKG)   , intent(in)    :: x
    !>                                                  integer(IK) , intent(in)    :: order
    !>                                                  real(RKG)                   :: func
    !>                                              end function
    !>                                          \endcode
    !>                                          where `RKG` refers to the kind of the output argument `root`.<br>
    !>                              </ol>
    !>  \param[inout]   root    :   The input/output scalar of type `real` of kind \RKALL.<br>
    !>                              On input,<br>
    !>                              <ol>
    !>                                  <li>    If the specified `method` is of type [newton_type](@ref pm_mathRoot::newton_type) or [halley_type](@ref pm_mathRoot::halley_type),
    !>                                          and if `lb < root .and. root <ub`, then the specified input value for `root` will be used as the best guess initial value
    !>                                          for the function `root` within the bracket specified by the input arguments `lb` and `ub`.<br>
    !>                                          Otherwise the center of the specifie dinput bracket `[lb, ub]` wil be used as the initial best-guess.<br>
    !>                                  <li>    Any input value for `root` is ignored for all other root-finding methods.<br>
    !>                              </ol>
    !>                              On output, `root` contains the best approximation to the root of the user-specified
    !>                              function `getFunc()` within the specified search interval `[lb, ub]` and tolerance threshold `abstol`.<br>
    !>  \param[in]      lb      :   The input scalar of same type and kind as the output `root`, representing the lower bound of the bracket (search) interval.<br>
    !>  \param[in]      ub      :   The input scalar of same type and kind as the output `root`, representing the upper bound of the bracket (search) interval.<br>
    !>  \param[in]      abstol  :   The input scalar of same type and kind as the output `root`, representing the absolute tolerance used as the stopping criterion of the search.<br>
    !>                              The iterations of the specified input `method` continue until the search interval becomes smaller than `abstol` in absolute units.<br>
    !>                              Care must be taken for specifying a reasonable value for `abstol` (see the warnings below).<br>
    !>                              If no suitable value for `abstol` is known a priori, try `abstol = epsilon(0._RKG)**.8 * (abs(lb) + abs(ub))`
    !>                              where `RKG` refers to the kind of the output argument `root`.<br>
    !>  \param[out]     neval   :   The output scalar argument of type `integer` of default kind \IK, containing the total number of `getFunc()` function calls made by the algorithm.<br>
    !>                              <ol>
    !>                                  <li>    A positive output `neval` implies the **successful** convergence of the algorithm to the function root after `+neval` function evaluations.<br>
    !>                                  <li>    A negative output `neval` implies the **failure** of convergence of the algorithm to the function root after `-neval` function evaluations.<br>
    !>                                  <li>    A zero output `neval` only occurs if the root lies at the search bracket boundaries specified by the input argument `lb` or `ub`.<br>
    !>                              </ol>
    !>                              Convergence failures rarely occur. If they do, setting the input arguments `abstol` and/or `niter` to larger values may resolve the failure.<br>
    !>  \param[in]      niter   :   The input scalar of type `integer` of default kind \IK, representing the **maximum number of iterations allowed** for the specified `method` to find the root.<br>
    !>                              The default number of steps `niter` within the algorithm is a compile-time constant that depends on the `kind` of the `real` input arguments.<br>
    !>                              (**optional**, default = `ceiling(2 * precision(lb) / log10(2.))`)
    !>
    !>  \interface{setRoot}
    !>  \code{.F90}
    !>
    !>      use pm_mathRoot, only: setRoot
    !>
    !>      call setRoot(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
    !>      call setRoot(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `lb < ub` must hold for the corresponding input arguments.<br>
    !>  The condition `0._RKG < abstol` must hold for the corresponding input arguments.<br>
    !>  The condition `abstol < (ub - lb)` must hold for the corresponding input arguments.<br>
    !>  The condition `abs(getFunc(lb), lf) < abstol` must hold for the corresponding input arguments.<br>
    !>  The condition `abs(getFunc(ub), uf) < abstol` must hold for the corresponding input arguments.<br>
    !>  The condition `sign(lf, 1.) /= sign(uf, 1.)` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < niter` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warning
    !>  It is crucial to keep in mind that computers use a fixed number of binary digits to represent floating-point numbers.<br>
    !>  While the user-specified function might analytically pass through zero, it is possible that its computed value is never zero, for any floating-point argument.<br>
    !>  One must also decide on what accuracy on the root is attainable.<br>
    !>  For example, convergence to within \f$10^{−6}\f$ in absolute value is reasonable when the root lies near \f$1\f$,
    !>  but unachievable if the root lies near an extremely large number near floating-point overflow.<br>
    !>
    !>  \impure
    !>
    !>  \recursive
    !>
    !>  \see
    !>  [getRoot](@ref pm_mathRoot::getRoot)<br>
    !>  [setRoot](@ref pm_mathRoot::setRoot)<br>
    !>  [getPolyRoot](@ref pm_polynomial::getPolyRoot)<br>
    !>  [setPolyRoot](@ref pm_polynomial::setPolyRoot)<br>
    !>
    !>  \example{setRoot}
    !>  \include{lineno} example/pm_mathRoot/setRoot/main.F90
    !>  \compilef{setRoot}
    !>  \output{setRoot}
    !>  \include{lineno} example/pm_mathRoot/setRoot/main.out.F90
    !>
    !>  \test
    !>  [test_pm_mathRoot](@ref test_pm_mathRoot)
    !>
    !>  \todo
    !>  \phigh
    !>  The current implementation of the Schroder method does not converge for **certain problems** that the Halley method converges readily.<br>
    !>  This may be due to the nature of the Schroder update.<br>
    !>  Regardless, this should be further investigated.<br>
    !>  For now, Schroder is diverted to Halley.<br>
    !>
    !>  \final{setRoot}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

    ! False

    interface setRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module subroutine setRootFalseFixed_RK5(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootFalseFixed_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(false_type)        , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK4_ENABLED
    recursive module subroutine setRootFalseFixed_RK4(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootFalseFixed_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(false_type)        , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK3_ENABLED
    recursive module subroutine setRootFalseFixed_RK3(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootFalseFixed_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(false_type)        , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK2_ENABLED
    recursive module subroutine setRootFalseFixed_RK2(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootFalseFixed_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(false_type)        , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK1_ENABLED
    recursive module subroutine setRootFalseFixed_RK1(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootFalseFixed_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(false_type)        , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module subroutine setRootFalseNiter_RK5(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootFalseNiter_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(false_type)        , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK4_ENABLED
    recursive module subroutine setRootFalseNiter_RK4(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootFalseNiter_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(false_type)        , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK3_ENABLED
    recursive module subroutine setRootFalseNiter_RK3(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootFalseNiter_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(false_type)        , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK2_ENABLED
    recursive module subroutine setRootFalseNiter_RK2(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootFalseNiter_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(false_type)        , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK1_ENABLED
    recursive module subroutine setRootFalseNiter_RK1(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootFalseNiter_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(false_type)        , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Bisection

    interface setRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module subroutine setRootBisectionFixed_RK5(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootBisectionFixed_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(bisection_type)    , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK4_ENABLED
    recursive module subroutine setRootBisectionFixed_RK4(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootBisectionFixed_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(bisection_type)    , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK3_ENABLED
    recursive module subroutine setRootBisectionFixed_RK3(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootBisectionFixed_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(bisection_type)    , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK2_ENABLED
    recursive module subroutine setRootBisectionFixed_RK2(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootBisectionFixed_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(bisection_type)    , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK1_ENABLED
    recursive module subroutine setRootBisectionFixed_RK1(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootBisectionFixed_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(bisection_type)    , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module subroutine setRootBisectionNiter_RK5(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootBisectionNiter_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(bisection_type)    , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK4_ENABLED
    recursive module subroutine setRootBisectionNiter_RK4(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootBisectionNiter_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(bisection_type)    , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK3_ENABLED
    recursive module subroutine setRootBisectionNiter_RK3(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootBisectionNiter_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(bisection_type)    , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK2_ENABLED
    recursive module subroutine setRootBisectionNiter_RK2(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootBisectionNiter_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(bisection_type)    , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK1_ENABLED
    recursive module subroutine setRootBisectionNiter_RK1(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootBisectionNiter_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(bisection_type)    , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Secant

    interface setRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module subroutine setRootSecantFixed_RK5(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootSecantFixed_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(secant_type)       , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK4_ENABLED
    recursive module subroutine setRootSecantFixed_RK4(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootSecantFixed_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(secant_type)       , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK3_ENABLED
    recursive module subroutine setRootSecantFixed_RK3(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootSecantFixed_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(secant_type)       , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK2_ENABLED
    recursive module subroutine setRootSecantFixed_RK2(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootSecantFixed_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(secant_type)       , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK1_ENABLED
    recursive module subroutine setRootSecantFixed_RK1(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootSecantFixed_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(secant_type)       , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module subroutine setRootSecantNiter_RK5(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootSecantNiter_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(secant_type)       , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK4_ENABLED
    recursive module subroutine setRootSecantNiter_RK4(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootSecantNiter_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(secant_type)       , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK3_ENABLED
    recursive module subroutine setRootSecantNiter_RK3(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootSecantNiter_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(secant_type)       , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK2_ENABLED
    recursive module subroutine setRootSecantNiter_RK2(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootSecantNiter_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(secant_type)       , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK1_ENABLED
    recursive module subroutine setRootSecantNiter_RK1(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootSecantNiter_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(secant_type)       , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Brent

    interface setRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module subroutine setRootBrentFixed_RK5(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootBrentFixed_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(brent_type)        , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK4_ENABLED
    recursive module subroutine setRootBrentFixed_RK4(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootBrentFixed_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(brent_type)        , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK3_ENABLED
    recursive module subroutine setRootBrentFixed_RK3(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootBrentFixed_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(brent_type)        , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK2_ENABLED
    recursive module subroutine setRootBrentFixed_RK2(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootBrentFixed_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(brent_type)        , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK1_ENABLED
    recursive module subroutine setRootBrentFixed_RK1(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootBrentFixed_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(brent_type)        , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module subroutine setRootBrentNiter_RK5(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootBrentNiter_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(brent_type)        , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK4_ENABLED
    recursive module subroutine setRootBrentNiter_RK4(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootBrentNiter_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(brent_type)        , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK3_ENABLED
    recursive module subroutine setRootBrentNiter_RK3(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootBrentNiter_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(brent_type)        , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK2_ENABLED
    recursive module subroutine setRootBrentNiter_RK2(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootBrentNiter_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(brent_type)        , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK1_ENABLED
    recursive module subroutine setRootBrentNiter_RK1(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootBrentNiter_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(brent_type)        , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Ridders

    interface setRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module subroutine setRootRiddersFixed_RK5(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootRiddersFixed_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(ridders_type)      , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK4_ENABLED
    recursive module subroutine setRootRiddersFixed_RK4(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootRiddersFixed_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(ridders_type)      , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK3_ENABLED
    recursive module subroutine setRootRiddersFixed_RK3(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootRiddersFixed_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(ridders_type)      , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK2_ENABLED
    recursive module subroutine setRootRiddersFixed_RK2(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootRiddersFixed_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(ridders_type)      , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK1_ENABLED
    recursive module subroutine setRootRiddersFixed_RK1(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootRiddersFixed_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(ridders_type)      , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module subroutine setRootRiddersNiter_RK5(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootRiddersNiter_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(ridders_type)      , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK4_ENABLED
    recursive module subroutine setRootRiddersNiter_RK4(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootRiddersNiter_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(ridders_type)      , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK3_ENABLED
    recursive module subroutine setRootRiddersNiter_RK3(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootRiddersNiter_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(ridders_type)      , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK2_ENABLED
    recursive module subroutine setRootRiddersNiter_RK2(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootRiddersNiter_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(ridders_type)      , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK1_ENABLED
    recursive module subroutine setRootRiddersNiter_RK1(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootRiddersNiter_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(ridders_type)      , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! TOMS748

    interface setRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module subroutine setRootTOMS748Fixed_RK5(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootTOMS748Fixed_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(toms748_type)      , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK4_ENABLED
    recursive module subroutine setRootTOMS748Fixed_RK4(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootTOMS748Fixed_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(toms748_type)      , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK3_ENABLED
    recursive module subroutine setRootTOMS748Fixed_RK3(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootTOMS748Fixed_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(toms748_type)      , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK2_ENABLED
    recursive module subroutine setRootTOMS748Fixed_RK2(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootTOMS748Fixed_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(toms748_type)      , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK1_ENABLED
    recursive module subroutine setRootTOMS748Fixed_RK1(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootTOMS748Fixed_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(toms748_type)      , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module subroutine setRootTOMS748Niter_RK5(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootTOMS748Niter_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(toms748_type)      , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK4_ENABLED
    recursive module subroutine setRootTOMS748Niter_RK4(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootTOMS748Niter_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(toms748_type)      , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK3_ENABLED
    recursive module subroutine setRootTOMS748Niter_RK3(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootTOMS748Niter_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(toms748_type)      , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK2_ENABLED
    recursive module subroutine setRootTOMS748Niter_RK2(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootTOMS748Niter_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(toms748_type)      , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK1_ENABLED
    recursive module subroutine setRootTOMS748Niter_RK1(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootTOMS748Niter_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(toms748_type)      , intent(in)                :: method
        real(RKG)               , intent(out)               :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Newton

    interface setRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module subroutine setRootNewtonFixed_RK5(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootNewtonFixed_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(newton_type)       , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK4_ENABLED
    recursive module subroutine setRootNewtonFixed_RK4(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootNewtonFixed_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(newton_type)       , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK3_ENABLED
    recursive module subroutine setRootNewtonFixed_RK3(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootNewtonFixed_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(newton_type)       , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK2_ENABLED
    recursive module subroutine setRootNewtonFixed_RK2(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootNewtonFixed_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(newton_type)       , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK1_ENABLED
    recursive module subroutine setRootNewtonFixed_RK1(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootNewtonFixed_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(newton_type)       , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module subroutine setRootNewtonNiter_RK5(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootNewtonNiter_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(newton_type)       , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK4_ENABLED
    recursive module subroutine setRootNewtonNiter_RK4(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootNewtonNiter_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(newton_type)       , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK3_ENABLED
    recursive module subroutine setRootNewtonNiter_RK3(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootNewtonNiter_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(newton_type)       , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK2_ENABLED
    recursive module subroutine setRootNewtonNiter_RK2(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootNewtonNiter_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(newton_type)       , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK1_ENABLED
    recursive module subroutine setRootNewtonNiter_RK1(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootNewtonNiter_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(newton_type)       , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Halley

    interface setRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module subroutine setRootHalleyFixed_RK5(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootHalleyFixed_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(halley_type)       , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK4_ENABLED
    recursive module subroutine setRootHalleyFixed_RK4(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootHalleyFixed_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(halley_type)       , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK3_ENABLED
    recursive module subroutine setRootHalleyFixed_RK3(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootHalleyFixed_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(halley_type)       , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK2_ENABLED
    recursive module subroutine setRootHalleyFixed_RK2(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootHalleyFixed_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(halley_type)       , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK1_ENABLED
    recursive module subroutine setRootHalleyFixed_RK1(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootHalleyFixed_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(halley_type)       , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module subroutine setRootHalleyNiter_RK5(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootHalleyNiter_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(halley_type)       , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK4_ENABLED
    recursive module subroutine setRootHalleyNiter_RK4(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootHalleyNiter_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(halley_type)       , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK3_ENABLED
    recursive module subroutine setRootHalleyNiter_RK3(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootHalleyNiter_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(halley_type)       , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK2_ENABLED
    recursive module subroutine setRootHalleyNiter_RK2(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootHalleyNiter_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(halley_type)       , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK1_ENABLED
    recursive module subroutine setRootHalleyNiter_RK1(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootHalleyNiter_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(halley_type)       , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Schroder

    interface setRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module subroutine setRootSchroderFixed_RK5(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootSchroderFixed_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(schroder_type)     , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK4_ENABLED
    recursive module subroutine setRootSchroderFixed_RK4(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootSchroderFixed_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(schroder_type)     , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK3_ENABLED
    recursive module subroutine setRootSchroderFixed_RK3(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootSchroderFixed_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(schroder_type)     , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK2_ENABLED
    recursive module subroutine setRootSchroderFixed_RK2(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootSchroderFixed_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(schroder_type)     , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

#if RK1_ENABLED
    recursive module subroutine setRootSchroderFixed_RK1(method, getFunc, root, lb, ub, lf, uf, abstol, neval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootSchroderFixed_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(schroder_type)     , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module subroutine setRootSchroderNiter_RK5(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootSchroderNiter_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        type(schroder_type)     , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK4_ENABLED
    recursive module subroutine setRootSchroderNiter_RK4(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootSchroderNiter_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        type(schroder_type)     , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK3_ENABLED
    recursive module subroutine setRootSchroderNiter_RK3(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootSchroderNiter_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        type(schroder_type)     , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK2_ENABLED
    recursive module subroutine setRootSchroderNiter_RK2(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootSchroderNiter_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        type(schroder_type)     , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

#if RK1_ENABLED
    recursive module subroutine setRootSchroderNiter_RK1(method, getFunc, root, lb, ub, lf, uf, abstol, neval, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRootSchroderNiter_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        type(schroder_type)     , intent(in)                :: method
        real(RKG)               , intent(inout)             :: root
        real(RKG)               , value                     :: lb, ub
        real(RKG)               , value                     :: lf, uf
        real(RKG)               , intent(in)                :: abstol
        integer(IK)             , intent(out)               :: neval
        integer(IK)             , intent(in)                :: niter
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_mathRoot