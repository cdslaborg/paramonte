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
!>  This module contains procedures and generic interfaces for performing various mathematical operations involving polynomials.<br>
!>
!>  \details
!>  Specifically, this module contains generic interfaces for the following computational polynomial tasks:<br>
!>  <ol>
!>      <li>    Evaluation of [polynomial functions](@ref pm_polynomial::getPolyVal) using the Horner method.<br>
!>      <li>    Evaluation of the [roots of polynomials](@ref pm_polynomial::getPolyRoot) of arbitrary degrees.<br>
!>      <li>    Evaluation of the [addition of two polynomials](@ref pm_polynomial::getPolyAdd) of arbitrary degrees.<br>
!>      <li>    Evaluation of the [subtraction of two polynomials](@ref pm_polynomial::getPolySub) of arbitrary degrees.<br>
!>      <li>    Evaluation of the [differentiation of polynomials](@ref pm_polynomial::getPolyDiff) of arbitrary degrees.<br>
!>      <li>    converting univariate polynomial coefficients to [readable univariate polynomial expressions in string format](@ref pm_polynomial::getPolyStr) for better display.<br>
!>  </ol>
!>
!>  Introduction
!>  ------------
!>
!>  In mathematics, a polynomial is an expression consisting of **indeterminates** (also called **variables**) and **coefficients**,
!>  that involves only the operations of **addition**, **subtraction**, **multiplication**, and **positive-integer powers** of variables.<br>
!>  An example of a polynomial of a single indeterminate \f$x\f$ is \f$x^2 − 4x + 7\f$.<br>
!>  An example with three indeterminates is \f$x^3 + 2xyz2 − yz + 1\f$.<br>
!>  Polynomials appear in many areas of mathematics and science.<br>
!>  For example, they are used to form polynomial equations, which encode a wide range of problems.<br>
!>  They are used to define polynomial functions, which appear in settings ranging from basic chemistry and physics to economics and social science.<br>
!>  They are used in calculus and numerical analysis to approximate other functions.<br>
!>  In advanced mathematics, polynomials are used to construct polynomial rings and algebraic varieties, which are central concepts in algebra and algebraic geometry.<br>
!>
!>  Definition
!>  ----------
!>
!>  A polynomial expression is an expression that can be built from constants and symbols called variables or
!>  indeterminates by means of addition, multiplication and exponentiation to a non-negative integer power.<br>
!>  The constants are generally numbers, but may be any expression that do not involve the indeterminates,
!>  and represent mathematical objects that can be added and multiplied.<br>
!>  Two polynomial expressions are considered as defining the same polynomial if they may be transformed,
!>  one to the other, by applying the usual properties of commutativity,
!>  associativity and distributivity of addition and multiplication.<br>
!>  For example \f$(x-1)(x-2)\f$ and \f$x^{2}-3x+2\f$ are two polynomial expressions that represent the same polynomial.<br>
!>  Therefore, one has the equality \f$(x-1)(x-2) = x^{2}-3x+2\f$.<br>
!>  A polynomial in a single indeterminate x can always be written (or rewritten) in the form
!>  \f{equation}{
!>      a_{n}x^{n}+a_{n-1}x^{n-1}+\dotsb +a_{2}x^{2}+a_{1}x+a_{0} ~,
!>  \f}
!>  where \f$a_{0}, \ldots , a_{n}\f$ are constants that are called the coefficients of the polynomial, and \f$x\f$ is the indeterminate.<br>
!>  The word **indeterminate** means that \f$x\f$ represents no particular value, although any value may be substituted for it.<br>
!>  The mapping that associates the result of this substitution to the substituted value is a function, called a **polynomial function**.<br>
!>  This can be expressed more concisely by using summation notation:<br>
!>  \f{equation}{
!>      \sum_{k = 0}^{n} a_{k}x^{k} ~.
!>  \f}
!>
!>  In other words, a polynomial can either be zero or can be written as the sum of a finite number of non-zero terms.<br>
!>  Each term consists of the product of a number - called the coefficient of the term – and a finite number of indeterminates, raised to non-negative integer powers.<br>
!>
!>  Addition
!>  --------
!>
!>  Polynomials can be added using the associative law of addition (grouping all their terms together into a single sum),
!>  possibly followed by reordering (using the commutative law) and combining of like terms.<br>
!>  For example, if \f$P = 3x^{2} - 2x + 5xy - 2\f$ and \f$Q = -3x^{2} + 3x + 4y^{2} + 8\f$, then the sum,
!>  \f{equation}{
!>      P + Q = 3x^{2} - 2x + 5xy - 2 - 3x^{2} + 3x + 4y^{2} + 8 ~,
!>  \f}
!>  can be reordered and regrouped as,
!>  \f{equation}{
!>      P + Q = (3x^{2} - 3x^{2}) + (-2x + 3x) + 5xy + 4y^{2} + (8 - 2) ~,
!>  \f}
!>  and then simplified to,
!>  \f{equation}{
!>      P + Q = x + 5xy + 4y^{2} + 6 ~.
!>  \f}
!>
!>  When polynomials are added together, the result is another polynomial.<br>
!>  In summary, to add two polynomials, simply add the corresponding coefficients of terms of equal exponent in the two polynomials.<br>
!>  This results in a new set of polynomial coefficients which represent the polynomial resulting from the sum of the two polynomials.<br>
!>
!>  \note
!>  The degree of the resulting polynomial from the addition of two other
!>  polynomials is always the maximum of the degrees of the two polynomials added.<br>
!>
!>  Subtraction
!>  -----------
!>
!>  Subtraction of polynomials is similar to their additions.<br>
!>  To subtract two polynomials, simply subtract the corresponding coefficients of terms of equal exponent in the two polynomials.<br>
!>  This results in a new set of polynomial coefficients which represent the polynomial resulting from the subtraction of the two polynomials.<br>
!>
!>  \note
!>  The degree of the resulting polynomial from the subtraction of two other polynomials is always
!>  the maximum of the degrees of the two polynomials subtracted.<br>
!>
!>  Multiplication
!>  --------------
!>
!>  Polynomials can also be multiplied.<br>
!>  To expand the product of two polynomials into a sum of terms, the distributive law is repeatedly applied,
!>  which results in each term of one polynomial being multiplied by every term of the other.<br>
!>  For example, if
!>  \f{equation}{
!>      \begin{aligned}
!>      \color{Red}     P & \color{Red}     {= 2x + 3y + 5} \\
!>      \color{Blue}    Q & \color{Blue}    {= 2x + 5y + xy + 1}
!>  \end{aligned}
!>  \f}
!>  then,
!>  \f{equation}{
!>      \begin{array}{rccrcrcrcr}
!>          {\color {Red}{P}}{\color {Blue}{Q}}&{=}&&({\color {Red}{2x}}\cdot {\color {Blue}{2x}})& + &({\color {Red}{2x}}\cdot {\color {Blue}{5y}})& + &({\color {Red}{2x}}\cdot {\color {Blue}{xy}})& + &({\color {Red}{2x}}\cdot {\color {Blue}{1}}) \\
!>          && + &({\color {Red}{3y}}\cdot {\color {Blue}{2x}})&+&({\color {Red}{3y}}\cdot {\color {Blue}{5y}})&+&({\color {Red}{3y}}\cdot {\color {Blue}{xy}})&+&({\color {Red}{3y}}\cdot {\color {Blue}{1}}) \\
!>          && + &({\color {Red}{5}}\cdot {\color {Blue}{2x}})&+&({\color {Red}{5}}\cdot {\color {Blue}{5y}})&+&({\color {Red}{5}}\cdot {\color {Blue}{xy}})&+&({\color {Red}{5}}\cdot {\color {Blue}{1}})
!>      \end{array}
!>  \f}
!>
!>  Carrying out the multiplication in each term produces,
!>  \f{equation}{
!>      \begin{array}{rccrcrcrcr}
!>          PQ &=&& 4x^{2}&+&10xy&+&2x^{2}y&+&2x \\
!>          &&+&6xy&+&15y^{2}&+&3xy^{2}&+&3y\\
!>          &&+&10x&+&25y&+&5xy&+&5.
!>      \end{array}
!>  \f}
!>
!>  Combining similar terms yields,
!>  \f{equation}{
!>      \begin{array}{rcccrcrcrcr}
!>          PQ &=&& 4x^{2}&+&(10xy+6xy+5xy)&+&2x^{2}y&+&(2x+10x) \\
!>          &&+&15y^{2}&+&3xy^{2}&+&(3y+25y)&+&5
!>      \end{array}
!>  \f}
!>  which can be simplified to,
!>  \f{equation}{
!>  PQ = 4x^{2}+21xy+2x^{2}y+12x+15y^{2}+3xy^{2}+28y+5 ~.
!>  \f}
!>
!>  As in the example, **the product of polynomials is always a polynomial**.<br>
!>
!>  In summary,
!>  To expand the product of two polynomials into a sum of terms, the distributive law is repeatedly applied.<br>
!>  This results in each term of one polynomial being multiplied by every term of the other.<br>
!>
!>  \note
!>  The degree of the resulting polynomial from the multiplication of two other polynomials is always
!>  the sum of the degrees of the two multiplicands plus one.<br>
!>
!>  Division
!>  --------
!>
!>  The division of one polynomial by another **is not** typically a polynomial.<br>
!>  Instead, such ratios are a more general family of objects, called rational fractions, rational expressions, or rational functions, depending on context.<br>
!>  This is analogous to the fact that the ratio of two integers is a rational number, not necessarily an integer.<br>
!>  For example, the fraction \f$1 / (x^2 + 1)\f$ is not a polynomial, and it cannot be written as a finite sum of powers of the variable \f$x\f$.<br>
!>  For polynomials in one variable, there is a notion of **Euclidean division of polynomials**, generalizing the **Euclidean division of integers**.<br>
!>  This notion of the division \f$a(x) / b(x)\f$ results in two polynomials, a **quotient** \f$q(x)\f$ and a **remainder** \f$r(x)\f$, such that \f$a = b q + r\f$ and \f$\ms{degree}(r) < \ms{degree}(b)\f$.<br>
!>  The quotient and remainder may be computed by any of several algorithms, including **polynomial long division** and **synthetic division**.<br>
!>  When the denominator \f$b(x)\f$ is monic and linear, that is, \f$b(x) = x − c\f$ for some constant \f$c\f$,
!>  then the polynomial remainder theorem asserts that the remainder of the division of \f$a(x)\f$ by \f$b(x)\f$ is the evaluation \f$a(c)\f$.<br>
!>  In this case, the quotient may be computed by the Ruffini rule, a special case of synthetic division.<br>
!>
!>  Division algorithm
!>  ------------------
!>
!>  The methodology employed in this module relies on **polynomial long division**.<br>
!>  Polynomial long division is an algorithm for dividing a polynomial by another polynomial of the same or lower degree.<br>
!>  It is a generalized version of the familiar arithmetic technique called **long division**.<br>
!>  Another abbreviated method is **polynomial short division** (**Blomqvist method**).<br>
!>  The method of polynomial long division implements the Euclidean division of polynomials.<br>
!>  Starting from the polynomial \f$\ms{dividend}\f$ and a non-zero polynomial \f$\ms{divisor}\f$,
!>  it outputs a \f$\ms{Quotient}\f$ polynomial and a \f$\ms{Remainder}\f$ polynomial such that,
!>  \f{equation}{
!>      \ms{dividend} = \ms{divisor} \times \ms{Quotient} + \ms{Remainder} ~.
!>  \f}
!>  By definition, the resulting degree of \f$\ms{Remainder}\f$ is lower than the degree of \f$\ms{divisor}\f$.<br>
!>  The resulting polynomials \f$\ms{Quotient}\f$ and \f$\ms{Remainder}\f$ are unique and do not depend on the derivation methodology.<br>
!>  The result \f$\ms{Remainder} = 0\f$ occurs **if and only if** the polynomial \f$\ms{dividend}\f$ has \f$\ms{divisor}\f$ as a factor.<br>
!>
!>  \note
!>  The polynomial division is a means for testing whether one polynomial has another as a factor, and, if it does, for factoring it out.<br>
!>  For example, if \f$r\f$ is a root of the polynomial \f$\ms{dividend}\f$ is known, it can be factored out by dividing \f$\ms{dividend}\f$ by \f$\ms{divisor} = (x – r)\f$.<br>
!>
!>  [Polynomial long division](https://en.wikipedia.org/wiki/Polynomial_long_division)<br>
!>
!>  Derivative
!>  ----------
!>
!>  Calculating derivatives and integrals of polynomials is particularly simple, compared to other kinds of functions.<br>
!>  The derivative of the polynomial,
!>  \f{equation}{
!>      P = a_{n}x^{n}+a_{n-1}x^{n-1}+\dots +a_{2}x^{2}+a_{1}x+a_{0}=\sum _{i=0}^{n}a_{i}x^{i} ~,
!>  \f}
!>  with respect to \f$x\f$ is the polynomial,
!>  \f{equation}{
!>      na_{n}x^{n-1}+(n-1)a_{n-1}x^{n-2}+\dots +2a_{2}x+a_{1}=\sum _{i=1}^{n}ia_{i}x^{i-1} ~.
!>  \f}
!>
!>  Succinctly, given a polynomial of degree \f$\ms{degree}\f$,<br>
!>  \f{equation}{
!>      P(x) = \sum_{i ~=~ 0}^{i ~=~ \ms{degree}} ~ c_i x^i ~,
!>  \f}
!>  where \f$c_i\f$ is the \f$i\f$th coefficient of the polynomial, the \f$k\f$th-order derivative of the polynomial can be computed as,
!>  \f{equation}{
!>      \frac{\mathrm{d}^k P(x)}{\mathrm{d}x^k} =
!>      \begin{cases}
!>          \sum_{i ~=~ k}^{i ~=~ \ms{degree}} ~ \big( \prod_{j = i - k + 1}^{j = i} \big) c_i x^{(i - k)} ~,~& k < \ms{degree} \\
!>          0  ~,~& k \geq \ms{degree}
!>      \end{cases}
!>  \f}
!>
!>  \note
!>  The degree of the resulting polynomial from the \f$k\f$th-order differentiation of a
!>  polynomial of degree \f$\ms{degree}\f$ is by definition \f$\max(0, \ms{degree} - k)\f$.<br>
!>
!>  Integral
!>  --------
!>
!>  Similarly, the general **antiderivative** (or **indefinite integral**) of \f$P\f$ is,
!>  \f{equation}{
!>      {\frac {a_{n}x^{n+1}}{n+1}}+{\frac {a_{n-1}x^{n}}{n}}+\dots +{\frac {a_{2}x^{3}}{3}}+{\frac {a_{1}x^{2}}{2}}+a_{0}x+c=c+\sum _{i=0}^{n}{\frac {a_{i}x^{i+1}}{i+1}} ~,
!>  \f}
!>  where \f$c\f$ is an arbitrary constant.<br>
!>  For example, antiderivatives of \f$x^2 + 1\f$ have the form \f$\frac{1}{3}x^3 + x + c\f$.<br>
!>
!>  Polynomial functions and root-finding
!>  -------------------------------------
!>
!>  A polynomial equation, also called an algebraic equation, is an equation of the form,
!>  \f{equation}{
!>      a_{n}x^{n}+a_{n-1}x^{n-1}+\dotsb +a_{2}x^{2}+a_{1}x+a_{0}=0 ~.
!>  \f}
!>
!>  For example,<br>
!>  \f{equation}{
!>      3x^{2}+4x-5 = 0 ~,
!>  \f}
!>  is a polynomial equation.<br>
!>
!>  A **root** of a nonzero univariate polynomial \f$P\f$ is a value a of \f$x\f$ such that \f$P(a) = 0\f$.<br>
!>  In other words, a root of \f$P\f$ is a solution of the polynomial equation \f$P(x) = 0\f$ or a zero of the polynomial function defined by \f$P\f$.<br>
!>  In the case of the zero polynomial, every number is a zero of the corresponding function, and the concept of root is rarely considered.<br>
!>  A number \f$a\f$ is a root of a polynomial \f$P\f$ if and only if the linear polynomial \f$x − a\f$ divides \f$P\f$, that is if there is another polynomial \f$Q\f$ such that \f$P = (x − a) Q\f$.<br>
!>  It may happen that a power (greater than 1) of \f$x − a\f$ divides \f$P\f$; in this case, \f$a\f$ is a **multiple root** of \f$P\f$, and otherwise \f$a\f$ is a **simple root** of \f$P\f$.<br>
!>  If \f$P\f$ is a nonzero polynomial, there is a highest power \f$m\f$ such that \f$(x − a)m\f$ divides \f$P\f$, which is called the multiplicity of \f$a\f$ as a root of \f$P\f$.<br>
!>  The number of roots of a nonzero polynomial \f$P\f$, counted with their respective **multiplicities**, cannot exceed the degree of \f$P\f$,
!>  and equals this degree if all complex roots are considered (this is a consequence of the fundamental theorem of algebra).<br>
!>  The coefficients of a polynomial and its roots are related by the Vieta formulas.<br>
!>
!>  Some polynomials, such as \f$x^2 + 1\f$, do not have any roots among the real numbers.<br>
!>  If, however, the set of accepted solutions is expanded to the complex numbers, every non-constant polynomial has at least one root;<br>
!>  This is the fundamental theorem of algebra.<br>
!>  By successively dividing out factors \f$x − a\f$, one sees that any polynomial with complex coefficients can be written as a constant (its leading coefficient)
!>  times a product of such polynomial factors of degree \f$1\f$;<br>
!>  Consequently, **the number of (complex) roots counted with their multiplicities is exactly equal to the degree of a polynomial**.<br>
!>
!>  Polynomial root-finding: Eigenvalue method
!>  ------------------------------------------
!>
!>  Given the `real` or `complex` coefficients \f$c_i\f$ of a polynomial of degree \f$n > 0\f$,<br>
!>  \f{equation}{
!>      f(x) = \sum_{i = 0}^{n} ~ c_i x^{i} ~,
!>  \f}
!>  the procedures of this module compute the \f$n\f$ (potentially `complex`) roots of the polynomial \f$f(x)\f$.<br>
!>  The roots are always returned as `complex` numbers as polynomials with `real` coefficients can also have `complex` roots.<br>
!>  The algorithms of this module rely on finding the eigenvalues of a matrix \f$A\f$ that are the roots of its **characteristic polynomial** \f$f(x) = \mathrm{det}[A − xI]\f$.<br>
!>  It can be shown that a monic polynomial \f$f(x)\f$ is the characteristic polynomial of the
!>  special \f$n\times n\f$ square **Frobenius Companion Matrix**,<br>
!>  \f{equation}{
!>      C(f) =
!>      \begin{bmatrix}
!>          0 & 0 & \dots & 0 & -c_{0} \\
!>          1 & 0 & \dots & 0 & -c_{1} \\
!>          0 & 1 & \dots & 0 & -c_{2} \\
!>          \vdots & \vdots & \ddots & \vdots & \vdots \\
!>          0 & 0 & \dots & 1 & -c_{n-1}
!>      \end{bmatrix}
!>  \f}
!>  Although the eigenvalue root-finding method can be slower than other polynomial root-finding methods,
!>  it is generally a more robust technique than others.<br>
!>
!>  \see
!>  [pm_polynomial](@ref pm_polynomial)<br>
!>  [pm_polynomial](@ref pm_polynomial)<br>
!>  [MATH77](http://netlib.org/math)<br>
!>  [Root-Finding Algorithms](https://en.wikipedia.org/wiki/Root-finding_algorithms)<br>
!>  S. Goedecker, Remark on algorithms to find roots of polynomials, SIAM J. on Scientific Computing 15, 5 (Sept. 1994) 1058-1063.<br>
!>  B. S. Garbow, J. M. Boyle, J. J. Dongarra, and C. B. Moler, Matrix Eigensystem Routines | EISPACK Guide Extension, Lecture Notes in Computer Science 51, Springer Verlag, Berlin (1977) 343 pages.<br>
!>  B. T. Smith, J. M. Boyle, B. S. Garbow, Y. Ikebe, V. C. Klema, and C. B. Moler, Matrix Eigensystem Routines | EISPACK Guide, Lecture Notes in Computer Science 6, Springer Verlag, Berlin (1974) 387 pages.<br>
!>
!>  Polynomial root-finding: Jenkins–Traub method
!>  ---------------------------------------------
!>
!>  Given the `real` or `complex` coefficients \f$c_i\f$ of a polynomial of degree \f$n > 0\f$,<br>
!>  \f{equation}{
!>      f(x) = \sum_{i = 0}^{n} ~ c_i x^{i} ~,
!>  \f}
!>  the root-finding procedures of this module also compute the \f$n\f$ (potentially `complex`) roots of the polynomial \f$f(x)\f$ using the Jenkins–Traub method.<br>
!>  The roots are always returned as `complex` numbers as polynomials with `real` coefficients can also have `complex` roots.<br>
!>  The algorithms of this module rely on a three-stage algorithm for finding roots of polynomials with `real` or `complex` coefficients as outlined in Jenkins and Traub (1970).<br>
!>
!>  \note
!>  Although the Jenkins–Traub method is a popular fairly robust method of polynomial root-finding,
!>  it can still fail to find roots of certain highly unstable polynomials, some of which are exemplified in the original paper of Jenkins and Traub (1970).<br>
!>  In general, if the polynomial coefficients have highly variable number of significant digits, that is, they are of vastly different magnitudes (some very small and some very large),
!>  the Jenkins–Traub method might be prone to failure.<br>
!>  This can happen when, for example, the difference in the number of significant digits of smallest and largest coefficients is
!>  comparable to the precision of the `real` or `complex` kind used to represent the coefficients.<br>
!>  In such a case, using a higher-precision `real` or `complex` kind may resolve the instability.<br>
!>
!>  **Which polynomial root-finding method should I use?**<br>
!>  If you have a polynomial of highly varying coefficients, then [Eigenvalue method](@ref pm_polynomial) is likely going to be the more reliable approach.<br>
!>  The [Jenkins-Traub](@ref pm_polynomial) is also considered a relatively reliable fast **Sure-Fire** technique for finding the roots of polynomials.<br>
!>
!>  \remark
!>  The algorithms of this module that take polynomial coefficients of type `real` are expected to be
!>  generally up to four times faster than the corresponding algorithms that accept coefficients of type `complex`.<br>
!>  The extent of the accuracy of above claim is yet to be tested for the specific implementations of this module.<br>
!>
!>  \see
!>  [pm_mathRoot](@ref pm_mathRoot)<br>
!>  [Root-Finding Algorithms](https://en.wikipedia.org/wiki/Root-finding_algorithms)<br>
!>  Jenkins and Traub, 1970, *A Three-Stage Algorithm for Real Polynomials Using Quadratic Iteration*, SIAM.<br>
!>  Jenkins, ALGORITHM 493 - Zeros of a Real Polynomial [C2], ACM Transactions on Mathematical Software, Vol 1, No 2, June 1975.<br>
!>
!>  Polynomial root-finding: Laguerre method
!>  ----------------------------------------
!>
!>  In numerical analysis, the Laguerre method is a root-finding algorithm tailored to polynomials.<br>
!>  In other words, the Laguerre method can be used to numerically solve the equation \f$p(x) = 0\f$ for a given polynomial \f$p(x)\f$.<br>
!>  One of the most useful properties of this method is that it is, from extensive empirical study, very close to being a **sure-fire** method,
!>  meaning that it is **almost guaranteed to always converge to some root of the polynomial, no matter what initial guess is chosen**.<br>
!>  However, there are more efficient methods with which it is guaranteed to find all roots or all real roots.<br>
!>  This method is named in honor of the French mathematician [Edmond Laguerre](https://en.wikipedia.org/wiki/Edmond_Laguerre).<br>
!>
!>  The algorithm of the Laguerre method to find one root of a polynomial \f$p(x)\f$ of degree \f$n\f$ is:<br>
!>  <ol>
!>      <li>    Choose an initial guess \f$x_0\f$.
!>      <li>    For \f$k = 0, 1, 2, \ldots\f$,<br>
!>              <ol>
!>                  <li>    If \f$p(x_{k})\f$ is very small, exit the loop.<br>
!>                  <li>    Calculate \f$G = {\frac {p'(x_{k})}{p(x_{k})}}\f$.<br>
!>                  <li>    Calculate \f$H = G^{2} - {\frac {p''(x_{k})}{p(x_{k})}}\f$.<br>
!>                  <li>    Calculate \f$a = {\frac{n}{G\pm{\sqrt{(n - 1)(nH - G^{2})}}}}\f$, where
!>                          the sign is chosen to give the denominator with the larger absolute value, to avoid catastrophic cancellation.<br>
!>                  <li>    Set \f$x_{{k+1}}=x_{k} - a\f$.<br>
!>              </ol>
!>      <li>    Repeat until a is small enough or if the maximum number of iterations has been reached.
!>  </ol>
!>  If a root has been found, the corresponding linear factor can be removed from \f$p\f$.<br>
!>  This **deflation step** reduces the degree of the polynomial by one, so that eventually, approximations for all roots of \f$p\f$ can be found.<br>
!>  Note however that **deflation can lead to approximate factors that differ significantly from the corresponding exact factors**.<br>
!>  This error is least if the roots are found in the order of increasing magnitude.<br>
!>
!>  \benchmarks
!>
!>  \benchmark{setPolyRoot, The effects of `method` on runtime efficiency}
!>  The following program compares the runtime performance of [setPolyRoot](@ref pm_polynomial::setPolyRoot)
!>  using different polynomial root finding algorithms.<br>
!>
!>  \include{lineno} benchmark/pm_polynomial/setPolyRoot/main.F90
!>  \compilefb{setPolyRoot}
!>  \postprocb{setPolyRoot}
!>  \include{lineno} benchmark/pm_polynomial/setPolyRoot/main.py
!>  \visb{setPolyRoot}
!>  \image html benchmark/pm_polynomial/setPolyRoot/benchmark.setPolyRoot.runtime.png width=1000
!>  \image html benchmark/pm_polynomial/setPolyRoot/benchmark.setPolyRoot.runtime.ratio.png width=1000
!>  \image html benchmark/pm_polynomial/setPolyRoot/benchmark.setPolyRoot.root.count.png width=1000
!>  \moralb{setPolyRoot}
!>      -#  Among all root finding algorithms, [jenkins_type](@ref pm_polynomial::jenkins_type) appears to be the fastest.<br>
!>      -#  The [eigen_type](@ref pm_polynomial::eigen_type) method also tends to offer a comparably good performance.<br>
!>      -#  Unlike the above two, [laguerre_type](@ref pm_polynomial::laguerre_type) algorithm tends to significantly
!>          trail behind both in performance and reliability in finding all roots of the polynomial.<br>
!>
!>  \test
!>  [test_pm_polynomial](@ref test_pm_polynomial)<br>
!>
!>  \todo
!>  \phigh
!>  A generic interface `setPolyCoef(coef, root)` must be added that computes polynomial
!>  coefficients from its roots using recursive FFT-based polynomial multiplications.<br>
!>  See the commented-out generic interface `setPolyCoef` within this module as the starting point.<br>
!>
!>  \final
!>  Beware that the Skowron-Gould method of polynomial root finding in this module
!>  have the following LICENSE information carried over with them.<br>
!>
!>  \verbatim
!>      Copyright 2012 Jan Skowron & Andrew Gould
!>
!>      Licensed under the Apache License, Version 2.0 (the "License");
!>      you may not use this file except in compliance with the License.
!>      You may obtain a copy of the License at
!>
!>          http://www.apache.org/licenses/LICENSE-2.0
!>
!>      Unless required by applicable law or agreed to in writing, software
!>      distributed under the License is distributed on an "AS IS" BASIS,
!>      WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!>      See the License for the specific language governing permissions and
!>      limitations under the License.
!>
!>      -------------------------------------------------------------------!
!>
!>      The authors also make this file available under the terms of
!>      GNU Lesser General Public License version 2 or any later version.
!>      (text of the LGPL licence 2 in NOTICE file)
!>
!>      -------------------------------------------------------------------!
!>
!>      A custom in the scientific comunity is (regardless of the licence
!>      you chose to use or distribute this software under)
!>      that if this code was important in the scientific process or
!>      for the results of your scientific work, we kindly ask you for the
!>      appropriate citation of the Paper (Skowron & Gould 2012), and
!>      we would be greatful if you pass the information about
!>      the proper citation to anyone whom you redistribute this software to.
!>  \endverbatim
!>
!>  \author
!>  \FatemehBagheri, Tuesday 08:49 PM, August 10, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_polynomial

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_polynomial"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the value of the polynomial of arbitrary degree
    !>  whose coefficients are specified by the user in the order of increasing power.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_polynomial](@ref pm_polynomial) for details of the implementation.<br>
    !>
    !>  \param[in]  coef    :   The input array of,<br>
    !>                          <ul>
    !>                              <li>    type `complex` of the same kind as the kind of the input `x` of type `complex`, or<br>
    !>                              <li>    type `real` of the same kind as the kind of the input `x` of type `complex`, or<br>
    !>                              <li>    type `real` of the same kind as the kind of the input `x` of type `real`, <br>
    !>                          </ul>
    !>                          containing the coefficients of the polynomial in the order of increasing power.<br>
    !>                          By definition, the degree of the polynomial is `size(coef) - 1`.<br>
    !>  \param[in]  x       :   The input scalar of,<br>
    !>                          <ul>
    !>                              <li>    type `complex` of kind \CKALL,
    !>                              <li>    type `real` of kind \RKALL,
    !>                          </ul>
    !>                          representing the point at which the polynomial must be computed.
    !>
    !>  \return
    !>  `poly`              :   The output scalar of the same type and kind as the input `x`
    !>                          containing the polynomial value at the specified point.
    !>
    !>  \interface{getPolyVal}
    !>  \code{.F90}
    !>
    !>      use pm_polynomial, only: getPolyVal
    !>
    !>      poly = getPolyVal(coef(:), x)
    !>      poly(1:size(x)) = getPolyVal(coef(:), x(:))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < size(coef)` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \note
    !>  The current interface requirements of either,<br>
    !>  <ol>
    !>      <li>    `(complex :: coef, complex :: x)`, or<br>
    !>      <li>    `(real :: coef, complex :: x)`, or<br>
    !>      <li>    `(real :: coef, real :: x)`<br>
    !>  </ol>
    !>  are in line with the fact that polynomial coefficients of type `real` frequently have `complex` roots.<br>
    !>  In particular, the interface `(real :: coef, complex :: x)` allows seamless integration of the polynomial
    !>  evaluations of this module with various polynomial root-finding algorithm of the ParaMonte library.<br>
    !>
    !>  \see
    !>  [getPolyRoot](@ref pm_polynomial::getPolyRoot)<br>
    !>  [setPolyRoot](@ref pm_polynomial::setPolyRoot)<br>
    !>
    !>  \example{getPolyVal}
    !>  \include{lineno} example/pm_polynomial/getPolyVal/main.F90
    !>  \compilef{getPolyVal}
    !>  \output{getPolyVal}
    !>  \include{lineno} example/pm_polynomial/getPolyVal/main.out.F90
    !>
    !>  \test
    !>  [test_pm_polynomial](@ref test_pm_polynomial)
    !>
    !>  \final{getPolyVal}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 11:34 PM, August 10, 2021, Dallas, TX
    interface getPolyVal

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getPolyVal_D0_CK5_CK5(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D0_CK5_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG), intent(in)                    :: x
        complex(CKG)                                :: poly
    end function
#endif

#if CK4_ENABLED
    PURE module function getPolyVal_D0_CK4_CK4(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D0_CK4_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG), intent(in)                    :: x
        complex(CKG)                                :: poly
    end function
#endif

#if CK3_ENABLED
    PURE module function getPolyVal_D0_CK3_CK3(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D0_CK3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG), intent(in)                    :: x
        complex(CKG)                                :: poly
    end function
#endif

#if CK2_ENABLED
    PURE module function getPolyVal_D0_CK2_CK2(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D0_CK2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG), intent(in)                    :: x
        complex(CKG)                                :: poly
    end function
#endif

#if CK1_ENABLED
    PURE module function getPolyVal_D0_CK1_CK1(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D0_CK1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG), intent(in)                    :: x
        complex(CKG)                                :: poly
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getPolyVal_D0_RK5_CK5(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D0_RK5_CK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        complex(RKG), intent(in)                    :: x
        complex(RKG)                                :: poly
    end function
#endif

#if RK4_ENABLED
    PURE module function getPolyVal_D0_RK4_CK4(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D0_RK4_CK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        complex(RKG), intent(in)                    :: x
        complex(RKG)                                :: poly
    end function
#endif

#if RK3_ENABLED
    PURE module function getPolyVal_D0_RK3_CK3(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D0_RK3_CK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        complex(RKG), intent(in)                    :: x
        complex(RKG)                                :: poly
    end function
#endif

#if RK2_ENABLED
    PURE module function getPolyVal_D0_RK2_CK2(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D0_RK2_CK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        complex(RKG), intent(in)                    :: x
        complex(RKG)                                :: poly
    end function
#endif

#if RK1_ENABLED
    PURE module function getPolyVal_D0_RK1_CK1(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D0_RK1_CK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        complex(RKG), intent(in)                    :: x
        complex(RKG)                                :: poly
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getPolyVal_D0_RK5_RK5(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D0_RK5_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)   , intent(in)                    :: x
        real(RKG)                                   :: poly
    end function
#endif

#if RK4_ENABLED
    PURE module function getPolyVal_D0_RK4_RK4(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D0_RK4_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)   , intent(in)                    :: x
        real(RKG)                                   :: poly
    end function
#endif

#if RK3_ENABLED
    PURE module function getPolyVal_D0_RK3_RK3(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D0_RK3_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)   , intent(in)                    :: x
        real(RKG)                                   :: poly
    end function
#endif

#if RK2_ENABLED
    PURE module function getPolyVal_D0_RK2_RK2(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D0_RK2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)   , intent(in)                    :: x
        real(RKG)                                   :: poly
    end function
#endif

#if RK1_ENABLED
    PURE module function getPolyVal_D0_RK1_RK1(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D0_RK1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)   , intent(in)                    :: x
        real(RKG)                                   :: poly
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getPolyVal_D1_CK5_CK5(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D1_CK5_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG), intent(in)    , contiguous    :: x(:)
        complex(CKG)                                :: poly(size(x, 1, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getPolyVal_D1_CK4_CK4(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D1_CK4_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG), intent(in)    , contiguous    :: x(:)
        complex(CKG)                                :: poly(size(x, 1, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getPolyVal_D1_CK3_CK3(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D1_CK3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG), intent(in)    , contiguous    :: x(:)
        complex(CKG)                                :: poly(size(x, 1, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getPolyVal_D1_CK2_CK2(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D1_CK2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG), intent(in)    , contiguous    :: x(:)
        complex(CKG)                                :: poly(size(x, 1, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getPolyVal_D1_CK1_CK1(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D1_CK1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG), intent(in)    , contiguous    :: x(:)
        complex(CKG)                                :: poly(size(x, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getPolyVal_D1_RK5_CK5(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D1_RK5_CK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        complex(RKG), intent(in)    , contiguous    :: x(:)
        complex(RKG)                                :: poly(size(x, 1, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getPolyVal_D1_RK4_CK4(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D1_RK4_CK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        complex(RKG), intent(in)    , contiguous    :: x(:)
        complex(RKG)                                :: poly(size(x, 1, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getPolyVal_D1_RK3_CK3(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D1_RK3_CK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        complex(RKG), intent(in)    , contiguous    :: x(:)
        complex(RKG)                                :: poly(size(x, 1, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getPolyVal_D1_RK2_CK2(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D1_RK2_CK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        complex(RKG), intent(in)    , contiguous    :: x(:)
        complex(RKG)                                :: poly(size(x, 1, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getPolyVal_D1_RK1_CK1(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D1_RK1_CK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        complex(RKG), intent(in)    , contiguous    :: x(:)
        complex(RKG)                                :: poly(size(x, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getPolyVal_D1_RK5_RK5(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D1_RK5_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)   , intent(in)    , contiguous    :: x(:)
        real(RKG)                                   :: poly(size(x, 1, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getPolyVal_D1_RK4_RK4(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D1_RK4_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)   , intent(in)    , contiguous    :: x(:)
        real(RKG)                                   :: poly(size(x, 1, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getPolyVal_D1_RK3_RK3(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D1_RK3_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)   , intent(in)    , contiguous    :: x(:)
        real(RKG)                                   :: poly(size(x, 1, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getPolyVal_D1_RK2_RK2(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D1_RK2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)   , intent(in)    , contiguous    :: x(:)
        real(RKG)                                   :: poly(size(x, 1, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getPolyVal_D1_RK1_RK1(coef, x) result(poly)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyVal_D1_RK1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)   , intent(in)    , contiguous    :: x(:)
        real(RKG)                                   :: poly(size(x, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the vector of coefficients of the polynomial resulting from
    !>  the addition of a polynomial to another polynomial of arbitrary degrees.
    !>
    !>  \details
    !>  See the documentation of [pm_polynomial](@ref pm_polynomial) for details of the implementation.<br>
    !>
    !>  \param[in]  lhs     :   The input `contiguous` vector of non-zero size of,<br>
    !>                          <ul>
    !>                              <li>    type `complex` of kind \CKALL, or<br>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ul>
    !>                          containing the coefficients of the left-hand-side polynomial in the addition, in the order of increasing power.<br>
    !>                          By definition, the degree of the `lhs` polynomial is `size(lhs) - 1`.<br>
    !>                          This means that the condition `lhs(size(lhs)) /= 0.` must hold.<br>
    !>  \param[in]  rhs     :   The input `contiguous` vector of non-zero size of the same type and kind as `lhs`,
    !>                          containing the coefficients of the right-hand-side polynomial in the addition, in the order of increasing power.<br>
    !>                          By definition, the degree of the `rhs` polynomial is `size(rhs) - 1`.<br>
    !>                          This means that the condition `rhs(size(rhs)) /= 0.` must hold.<br>
    !>
    !>  \return
    !>  `add`               :   The output `contiguous` vector of the same type and kind as the input `lhs`, of size `max(size(lhs), size(rhs))`,
    !>                          containing the coefficients  (in the order of **increasing power**) of the resulting polynomial from the addition
    !>                          of the polynomial `lhs` of arbitrary degree to another polynomial `rhs` of arbitrary degree.<br>
    !>                          By definition, the degree of the `add` polynomial is `size(add) - 1`.<br>
    !>
    !>  \interface{getPolyAdd}
    !>  \code{.F90}
    !>
    !>      use pm_polynomial, only: getPolyAdd
    !>
    !>      add(1 : max(size(lhs), size(rhs))) = getPolyAdd(lhs(:), rhs(:))
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \note
    !>  An input empty `lhs` or `rhs` vector (of size zero) is interpreted as polynomial of degree zero with the single coefficient being `0`.<br>
    !>  This behavior is useful for flexibly inverting a polynomial division operation done via [setPolyDiv](@ref pm_polynomial::setPolyDiv).<br>
    !>
    !>  \see
    !>  [getPolyAdd](@ref pm_polynomial::getPolyAdd)<br>
    !>  [setPolyAdd](@ref pm_polynomial::setPolyAdd)<br>
    !>  [getPolySub](@ref pm_polynomial::getPolySub)<br>
    !>  [setPolySub](@ref pm_polynomial::setPolySub)<br>
    !>  [setPolyDiv](@ref pm_polynomial::setPolyDiv)<br>
    !>  [getPolyMul](@ref pm_polynomial::getPolyMul)<br>
    !>  [setPolyMul](@ref pm_polynomial::setPolyMul)<br>
    !>
    !>  \example{getPolyAdd}
    !>  \include{lineno} example/pm_polynomial/getPolyAdd/main.F90
    !>  \compilef{getPolyAdd}
    !>  \output{getPolyAdd}
    !>  \include{lineno} example/pm_polynomial/getPolyAdd/main.out.F90
    !>
    !>  \test
    !>  [test_pm_polynomial](@ref test_pm_polynomial)
    !>
    !>  \final{getPolyAdd}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 11:34 PM, August 10, 2021, Dallas, TX
    interface getPolyAdd

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getPolyAdd_CK5(lhs, rhs) result(add)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyAdd_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous    :: lhs(:)
        complex(CKG), intent(in)    , contiguous    :: rhs(:)
        complex(CKG)                                :: add(max(size(lhs), size(rhs)))
    end function
#endif

#if CK4_ENABLED
    PURE module function getPolyAdd_CK4(lhs, rhs) result(add)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyAdd_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: lhs(:)
        complex(CKG), intent(in)    , contiguous    :: rhs(:)
        complex(CKG)                                :: add(max(size(lhs), size(rhs)))
    end function
#endif

#if CK3_ENABLED
    PURE module function getPolyAdd_CK3(lhs, rhs) result(add)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyAdd_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: lhs(:)
        complex(CKG), intent(in)    , contiguous    :: rhs(:)
        complex(CKG)                                :: add(max(size(lhs), size(rhs)))
    end function
#endif

#if CK2_ENABLED
    PURE module function getPolyAdd_CK2(lhs, rhs) result(add)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyAdd_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: lhs(:)
        complex(CKG), intent(in)    , contiguous    :: rhs(:)
        complex(CKG)                                :: add(max(size(lhs), size(rhs)))
    end function
#endif

#if CK1_ENABLED
    PURE module function getPolyAdd_CK1(lhs, rhs) result(add)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyAdd_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: lhs(:)
        complex(CKG), intent(in)    , contiguous    :: rhs(:)
        complex(CKG)                                :: add(max(size(lhs), size(rhs)))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getPolyAdd_RK5(lhs, rhs) result(add)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyAdd_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: lhs(:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(:)
        real(RKG)                                   :: add(max(size(lhs), size(rhs)))
    end function
#endif

#if RK4_ENABLED
    PURE module function getPolyAdd_RK4(lhs, rhs) result(add)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyAdd_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: lhs(:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(:)
        real(RKG)                                   :: add(max(size(lhs), size(rhs)))
    end function
#endif

#if RK3_ENABLED
    PURE module function getPolyAdd_RK3(lhs, rhs) result(add)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyAdd_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: lhs(:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(:)
        real(RKG)                                   :: add(max(size(lhs), size(rhs)))
    end function
#endif

#if RK2_ENABLED
    PURE module function getPolyAdd_RK2(lhs, rhs) result(add)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyAdd_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: lhs(:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(:)
        real(RKG)                                   :: add(max(size(lhs), size(rhs)))
    end function
#endif

#if RK1_ENABLED
    PURE module function getPolyAdd_RK1(lhs, rhs) result(add)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyAdd_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: lhs(:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(:)
        real(RKG)                                   :: add(max(size(lhs), size(rhs)))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the vector of coefficients of the polynomial resulting from the addition of a polynomial to another polynomial of arbitrary degrees.
    !>
    !>  \details
    !>  See the documentation of [pm_polynomial](@ref pm_polynomial) for details of the implementation.<br>
    !>
    !>  \param[out] add     :   The output `contiguous` vector of the same type and kind as the input `lhs`, of size `max(size(lhs), size(rhs))`,
    !>                          containing the coefficients  (in the order of **increasing power**) of the resulting polynomial from the addition
    !>                          of the polynomial `lhs` of arbitrary degree to another polynomial `rhs` of arbitrary degree.<br>
    !>                          By definition, the degree of the `add` polynomial is `size(add) - 1`.<br>
    !>  \param[in]  lhs     :   The input `contiguous` vector of non-zero size of,<br>
    !>                          <ul>
    !>                              <li>    type `complex` of kind \CKALL, or<br>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ul>
    !>                          containing the coefficients of the left-hand-side polynomial in the addition, in the order of increasing power.<br>
    !>                          By definition, the degree of the `lhs` polynomial is `size(lhs) - 1`.<br>
    !>                          This means that the condition `lhs(size(lhs)) /= 0.` must hold.<br>
    !>  \param[in]  rhs     :   The input `contiguous` vector of non-zero size of the same type and kind as `lhs`,
    !>                          containing the coefficients of the right-hand-side polynomial in the addition, in the order of increasing power.<br>
    !>                          By definition, the degree of the `rhs` polynomial is `size(rhs) - 1`.<br>
    !>                          This means that the condition `rhs(size(rhs)) /= 0.` must hold.<br>
    !>
    !>  \interface{setPolyAdd}
    !>  \code{.F90}
    !>
    !>      use pm_polynomial, only: setPolyAdd
    !>
    !>      call setPolyAdd(add(1 : max(size(lhs), size(rhs))), lhs(:), rhs(:))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < size(rhs)` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < size(lhs)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(add) == max(size(lhs), size(rhs))` must hold for the corresponding input arguments.<br>
    !>  The condition `lhs(size(lhs)) /= 0.` must hold for the corresponding input arguments.<br>
    !>  The condition `rhs(size(rhs)) /= 0.` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getPolyAdd](@ref pm_polynomial::getPolyAdd)<br>
    !>  [setPolyAdd](@ref pm_polynomial::setPolyAdd)<br>
    !>  [getPolySub](@ref pm_polynomial::getPolySub)<br>
    !>  [setPolySub](@ref pm_polynomial::setPolySub)<br>
    !>  [setPolyDiv](@ref pm_polynomial::setPolyDiv)<br>
    !>  [getPolyMul](@ref pm_polynomial::getPolyMul)<br>
    !>  [setPolyMul](@ref pm_polynomial::setPolyMul)<br>
    !>
    !>  \example{setPolyAdd}
    !>  \include{lineno} example/pm_polynomial/setPolyAdd/main.F90
    !>  \compilef{setPolyAdd}
    !>  \output{setPolyAdd}
    !>  \include{lineno} example/pm_polynomial/setPolyAdd/main.out.F90
    !>
    !>  \test
    !>  [test_pm_polynomial](@ref test_pm_polynomial)
    !>
    !>  \final{setPolyAdd}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 11:34 PM, August 10, 2021, Dallas, TX
    interface setPolyAdd

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPolyAdd_CK5(add, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyAdd_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous    :: lhs(:)
        complex(CKG), intent(in)    , contiguous    :: rhs(:)
        complex(CKG), intent(out)   , contiguous    :: add(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPolyAdd_CK4(add, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyAdd_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: lhs(:)
        complex(CKG), intent(in)    , contiguous    :: rhs(:)
        complex(CKG), intent(out)   , contiguous    :: add(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPolyAdd_CK3(add, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyAdd_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: lhs(:)
        complex(CKG), intent(in)    , contiguous    :: rhs(:)
        complex(CKG), intent(out)   , contiguous    :: add(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPolyAdd_CK2(add, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyAdd_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: lhs(:)
        complex(CKG), intent(in)    , contiguous    :: rhs(:)
        complex(CKG), intent(out)   , contiguous    :: add(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPolyAdd_CK1(add, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyAdd_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: lhs(:)
        complex(CKG), intent(in)    , contiguous    :: rhs(:)
        complex(CKG), intent(out)   , contiguous    :: add(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPolyAdd_RK5(add, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyAdd_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: lhs(:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(:)
        real(RKG)   , intent(out)   , contiguous    :: add(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPolyAdd_RK4(add, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyAdd_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: lhs(:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(:)
        real(RKG)   , intent(out)   , contiguous    :: add(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPolyAdd_RK3(add, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyAdd_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: lhs(:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(:)
        real(RKG)   , intent(out)   , contiguous    :: add(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPolyAdd_RK2(add, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyAdd_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: lhs(:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(:)
        real(RKG)   , intent(out)   , contiguous    :: add(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPolyAdd_RK1(add, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyAdd_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: lhs(:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(:)
        real(RKG)   , intent(out)   , contiguous    :: add(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the vector of coefficients of the polynomial resulting from
    !>  the subtraction of a polynomial to another polynomial of arbitrary degrees.
    !>
    !>  \details
    !>  See the documentation of [pm_polynomial](@ref pm_polynomial) for details of the implementation.<br>
    !>
    !>  \param[in]  lhs     :   The input `contiguous` vector of non-zero size of,<br>
    !>                          <ul>
    !>                              <li>    type `complex` of kind \CKALL, or<br>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ul>
    !>                          containing the coefficients of the left-hand-side polynomial in the subtraction, in the order of increasing power.<br>
    !>                          By definition, the degree of the `lhs` polynomial is `size(lhs) - 1`.<br>
    !>                          This means that the condition `lhs(size(lhs)) /= 0.` must hold.<br>
    !>  \param[in]  rhs     :   The input `contiguous` vector of non-zero size of the same type and kind as `lhs`,
    !>                          containing the coefficients of the right-hand-side polynomial in the subtraction, in the order of increasing power.<br>
    !>                          By definition, the degree of the `rhs` polynomial is `size(rhs) - 1`.<br>
    !>                          This means that the condition `rhs(size(rhs)) /= 0.` must hold.<br>
    !>
    !>  \return
    !>  `sub`               :   The output `contiguous` vector of the same type and kind as the input `lhs`, of size `max(size(lhs), size(rhs))`,
    !>                          containing the coefficients  (in the order of **increasing power**) of the resulting polynomial from the subtraction
    !>                          of the polynomial `lhs` of arbitrary degree to another polynomial `rhs` of arbitrary degree.<br>
    !>                          By definition, the degree of the `sub` polynomial is `size(sub) - 1`.<br>
    !>
    !>  \interface{getPolySub}
    !>  \code{.F90}
    !>
    !>      use pm_polynomial, only: getPolySub
    !>
    !>      sub(1 : max(size(lhs), size(rhs))) = getPolySub(lhs(:), rhs(:))
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \note
    !>  An input empty `lhs` or `rhs` vector (of size zero) is interpreted as polynomial of degree zero with the single coefficient being `0`.<br>
    !>  This behavior is useful for flexibly inverting a polynomial division operation done via [setPolyDiv](@ref pm_polynomial::setPolyDiv).<br>
    !>
    !>  \see
    !>  [getPolySub](@ref pm_polynomial::getPolySub)<br>
    !>  [setPolySub](@ref pm_polynomial::setPolySub)<br>
    !>  [getPolySub](@ref pm_polynomial::getPolySub)<br>
    !>  [setPolySub](@ref pm_polynomial::setPolySub)<br>
    !>  [setPolyDiv](@ref pm_polynomial::setPolyDiv)<br>
    !>  [getPolyMul](@ref pm_polynomial::getPolyMul)<br>
    !>  [setPolyMul](@ref pm_polynomial::setPolyMul)<br>
    !>
    !>  \example{getPolySub}
    !>  \include{lineno} example/pm_polynomial/getPolySub/main.F90
    !>  \compilef{getPolySub}
    !>  \output{getPolySub}
    !>  \include{lineno} example/pm_polynomial/getPolySub/main.out.F90
    !>
    !>  \test
    !>  [test_pm_polynomial](@ref test_pm_polynomial)
    !>
    !>  \final{getPolySub}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 11:34 PM, August 10, 2021, Dallas, TX
    interface getPolySub

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getPolySub_CK5(lhs, rhs) result(sub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolySub_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous    :: lhs(:)
        complex(CKG), intent(in)    , contiguous    :: rhs(:)
        complex(CKG)                                :: sub(max(size(lhs), size(rhs)))
    end function
#endif

#if CK4_ENABLED
    PURE module function getPolySub_CK4(lhs, rhs) result(sub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolySub_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: lhs(:)
        complex(CKG), intent(in)    , contiguous    :: rhs(:)
        complex(CKG)                                :: sub(max(size(lhs), size(rhs)))
    end function
#endif

#if CK3_ENABLED
    PURE module function getPolySub_CK3(lhs, rhs) result(sub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolySub_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: lhs(:)
        complex(CKG), intent(in)    , contiguous    :: rhs(:)
        complex(CKG)                                :: sub(max(size(lhs), size(rhs)))
    end function
#endif

#if CK2_ENABLED
    PURE module function getPolySub_CK2(lhs, rhs) result(sub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolySub_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: lhs(:)
        complex(CKG), intent(in)    , contiguous    :: rhs(:)
        complex(CKG)                                :: sub(max(size(lhs), size(rhs)))
    end function
#endif

#if CK1_ENABLED
    PURE module function getPolySub_CK1(lhs, rhs) result(sub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolySub_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: lhs(:)
        complex(CKG), intent(in)    , contiguous    :: rhs(:)
        complex(CKG)                                :: sub(max(size(lhs), size(rhs)))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getPolySub_RK5(lhs, rhs) result(sub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolySub_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: lhs(:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(:)
        real(RKG)                                   :: sub(max(size(lhs), size(rhs)))
    end function
#endif

#if RK4_ENABLED
    PURE module function getPolySub_RK4(lhs, rhs) result(sub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolySub_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: lhs(:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(:)
        real(RKG)                                   :: sub(max(size(lhs), size(rhs)))
    end function
#endif

#if RK3_ENABLED
    PURE module function getPolySub_RK3(lhs, rhs) result(sub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolySub_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: lhs(:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(:)
        real(RKG)                                   :: sub(max(size(lhs), size(rhs)))
    end function
#endif

#if RK2_ENABLED
    PURE module function getPolySub_RK2(lhs, rhs) result(sub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolySub_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: lhs(:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(:)
        real(RKG)                                   :: sub(max(size(lhs), size(rhs)))
    end function
#endif

#if RK1_ENABLED
    PURE module function getPolySub_RK1(lhs, rhs) result(sub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolySub_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: lhs(:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(:)
        real(RKG)                                   :: sub(max(size(lhs), size(rhs)))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the vector of coefficients of the polynomial resulting from the subtraction of a polynomial to another polynomial of arbitrary degrees.
    !>
    !>  \details
    !>  See the documentation of [pm_polynomial](@ref pm_polynomial) for details of the implementation.<br>
    !>
    !>  \param[out] sub     :   The output `contiguous` vector of the same type and kind as the input `lhs`, of size `max(size(lhs), size(rhs))`,
    !>                          containing the coefficients  (in the order of **increasing power**) of the resulting polynomial from the subtraction
    !>                          of the polynomial `lhs` of arbitrary degree to another polynomial `rhs` of arbitrary degree.<br>
    !>                          By definition, the degree of the `sub` polynomial is `size(sub) - 1`.<br>
    !>  \param[in]  lhs     :   The input `contiguous` vector of non-zero size of,<br>
    !>                          <ul>
    !>                              <li>    type `complex` of kind \CKALL, or<br>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ul>
    !>                          containing the coefficients of the left-hand-side polynomial in the subtraction, in the order of increasing power.<br>
    !>                          By definition, the degree of the `lhs` polynomial is `size(lhs) - 1`.<br>
    !>                          This means that the condition `lhs(size(lhs)) /= 0.` must hold.<br>
    !>  \param[in]  rhs     :   The input `contiguous` vector of non-zero size of the same type and kind as `lhs`,
    !>                          containing the coefficients of the right-hand-side polynomial in the subtraction, in the order of increasing power.<br>
    !>                          By definition, the degree of the `rhs` polynomial is `size(rhs) - 1`.<br>
    !>                          This means that the condition `rhs(size(rhs)) /= 0.` must hold.<br>
    !>
    !>  \interface{setPolySub}
    !>  \code{.F90}
    !>
    !>      use pm_polynomial, only: setPolySub
    !>
    !>      call setPolySub(sub(1 : max(size(lhs), size(rhs))), lhs(:), rhs(:))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < size(rhs)` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < size(lhs)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(sub) == max(size(lhs), size(rhs))` must hold for the corresponding input arguments.<br>
    !>  The condition `lhs(size(lhs)) /= 0.` must hold for the corresponding input arguments.<br>
    !>  The condition `rhs(size(rhs)) /= 0.` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getPolySub](@ref pm_polynomial::getPolySub)<br>
    !>  [setPolySub](@ref pm_polynomial::setPolySub)<br>
    !>  [getPolySub](@ref pm_polynomial::getPolySub)<br>
    !>  [setPolySub](@ref pm_polynomial::setPolySub)<br>
    !>  [setPolyDiv](@ref pm_polynomial::setPolyDiv)<br>
    !>  [getPolyMul](@ref pm_polynomial::getPolyMul)<br>
    !>  [setPolyMul](@ref pm_polynomial::setPolyMul)<br>
    !>
    !>  \example{setPolySub}
    !>  \include{lineno} example/pm_polynomial/setPolySub/main.F90
    !>  \compilef{setPolySub}
    !>  \output{setPolySub}
    !>  \include{lineno} example/pm_polynomial/setPolySub/main.out.F90
    !>
    !>  \test
    !>  [test_pm_polynomial](@ref test_pm_polynomial)
    !>
    !>  \final{setPolySub}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 11:34 PM, August 10, 2021, Dallas, TX
    interface setPolySub

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPolySub_CK5(sub, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolySub_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous    :: lhs(:)
        complex(CKG), intent(in)    , contiguous    :: rhs(:)
        complex(CKG), intent(out)   , contiguous    :: sub(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPolySub_CK4(sub, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolySub_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: lhs(:)
        complex(CKG), intent(in)    , contiguous    :: rhs(:)
        complex(CKG), intent(out)   , contiguous    :: sub(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPolySub_CK3(sub, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolySub_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: lhs(:)
        complex(CKG), intent(in)    , contiguous    :: rhs(:)
        complex(CKG), intent(out)   , contiguous    :: sub(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPolySub_CK2(sub, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolySub_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: lhs(:)
        complex(CKG), intent(in)    , contiguous    :: rhs(:)
        complex(CKG), intent(out)   , contiguous    :: sub(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPolySub_CK1(sub, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolySub_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: lhs(:)
        complex(CKG), intent(in)    , contiguous    :: rhs(:)
        complex(CKG), intent(out)   , contiguous    :: sub(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPolySub_RK5(sub, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolySub_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: lhs(:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(:)
        real(RKG)   , intent(out)   , contiguous    :: sub(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPolySub_RK4(sub, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolySub_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: lhs(:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(:)
        real(RKG)   , intent(out)   , contiguous    :: sub(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPolySub_RK3(sub, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolySub_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: lhs(:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(:)
        real(RKG)   , intent(out)   , contiguous    :: sub(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPolySub_RK2(sub, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolySub_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: lhs(:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(:)
        real(RKG)   , intent(out)   , contiguous    :: sub(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPolySub_RK1(sub, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolySub_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: lhs(:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(:)
        real(RKG)   , intent(out)   , contiguous    :: sub(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the vector of coefficients of the polynomial resulting from
    !>  the multiplication of a polynomial with another polynomial of arbitrary degrees.
    !>
    !>  \details
    !>  See the documentation of [pm_polynomial](@ref pm_polynomial) for details of the implementation.<br>
    !>
    !>  \param[in]  lhs     :   The input `contiguous` vector of non-zero size of,<br>
    !>                          <ul>
    !>                              <li>    type `complex` of kind \CKALL, or<br>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ul>
    !>                          containing the coefficients of the left-hand-side polynomial in the multiplication, in the order of increasing power.<br>
    !>                          By definition, the degree of the `lhs` polynomial is `size(lhs) - 1`.<br>
    !>                          This means that the condition `lhs(size(lhs)) /= 0.` must hold.<br>
    !>  \param[in]  rhs     :   The input `contiguous` vector of non-zero size of the same type and kind as `lhs`,
    !>                          containing the coefficients of the right-hand-side polynomial in the multiplication, in the order of increasing power.<br>
    !>                          By definition, the degree of the `rhs` polynomial is `size(rhs) - 1`.<br>
    !>                          This means that the condition `rhs(size(rhs)) /= 0.` must hold.<br>
    !>
    !>  \return
    !>  `mul`               :   The output `contiguous` vector of the same type and kind as the input `lhs`, of size `max(0, min(size(lhs) * size(rhs), size(lhs) + size(rhs) - 1))`,
    !>                          containing the coefficients  (in the order of **increasing power**) of the resulting polynomial from the multiplication
    !>                          of the polynomial `lhs` of arbitrary degree with another polynomial `rhs` of arbitrary degree.<br>
    !>                          By definition, the degree of the `mul` polynomial is `size(mul) - 1`.<br>
    !>
    !>  \interface{getPolyMul}
    !>  \code{.F90}
    !>
    !>      use pm_polynomial, only: getPolyMul
    !>
    !>      mul(1 : size(lhs) + size(rhs) - 1) = getPolyMul(lhs(:), rhs(:))
    !>
    !>  \endcode
    !>
    !>  \note
    !>  The additional minimax criterion in size requirement of `mul(1:max(0, min(size(lhs) * size(rhs), size(lhs) + size(rhs) - 1)))`
    !>  is to ensure the that the procedures of this generic interface can gracefully handle input polynomials coefficient vectors `lhs` and `rhs` of zero sizes.<br>
    !>  A zero-sized polynomial coefficient vector is equivalent to the scalar zero.<br>
    !>  Such cases can occur in polynomial division operations where the output quotient or the remainder of the division is an empty vector.<br>
    !>  The specific size definition of `mul` enhances the flexibility and utility of the procedures of this generic interface in quick concise polynomial arithmetic.<br>
    !>  Note, however, that the same flexibility does not hold for the procedures of [setPolyMul](@ref pm_polynomial::setPolyMul) generic interface.<br>
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getPolySub](@ref pm_polynomial::getPolySub)<br>
    !>  [setPolySub](@ref pm_polynomial::setPolySub)<br>
    !>  [getPolySub](@ref pm_polynomial::getPolySub)<br>
    !>  [setPolySub](@ref pm_polynomial::setPolySub)<br>
    !>  [setPolyDiv](@ref pm_polynomial::setPolyDiv)<br>
    !>  [getPolyMul](@ref pm_polynomial::getPolyMul)<br>
    !>  [setPolyMul](@ref pm_polynomial::setPolyMul)<br>
    !>
    !>  \example{getPolyMul}
    !>  \include{lineno} example/pm_polynomial/getPolyMul/main.F90
    !>  \compilef{getPolyMul}
    !>  \output{getPolyMul}
    !>  \include{lineno} example/pm_polynomial/getPolyMul/main.out.F90
    !>
    !>  \test
    !>  [test_pm_polynomial](@ref test_pm_polynomial)
    !>
    !>  \final{getPolyMul}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 11:34 PM, August 10, 2021, Dallas, TX
    interface getPolyMul

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getPolyMul_CK5(lhs, rhs) result(mul)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyMul_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous    :: lhs(:)
        complex(CKG), intent(in)    , contiguous    :: rhs(:)
        complex(CKG)                                :: mul(1 : min(size(lhs, 1, IK) * size(rhs, 1, IK), size(lhs, 1, IK) + size(rhs, 1, IK) - 1_IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getPolyMul_CK4(lhs, rhs) result(mul)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyMul_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: lhs(:)
        complex(CKG), intent(in)    , contiguous    :: rhs(:)
        complex(CKG)                                :: mul(1 : min(size(lhs, 1, IK) * size(rhs, 1, IK), size(lhs, 1, IK) + size(rhs, 1, IK) - 1_IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getPolyMul_CK3(lhs, rhs) result(mul)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyMul_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: lhs(:)
        complex(CKG), intent(in)    , contiguous    :: rhs(:)
        complex(CKG)                                :: mul(1 : min(size(lhs, 1, IK) * size(rhs, 1, IK), size(lhs, 1, IK) + size(rhs, 1, IK) - 1_IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getPolyMul_CK2(lhs, rhs) result(mul)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyMul_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: lhs(:)
        complex(CKG), intent(in)    , contiguous    :: rhs(:)
        complex(CKG)                                :: mul(1 : min(size(lhs, 1, IK) * size(rhs, 1, IK), size(lhs, 1, IK) + size(rhs, 1, IK) - 1_IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getPolyMul_CK1(lhs, rhs) result(mul)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyMul_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: lhs(:)
        complex(CKG), intent(in)    , contiguous    :: rhs(:)
        complex(CKG)                                :: mul(1 : min(size(lhs, 1, IK) * size(rhs, 1, IK), size(lhs, 1, IK) + size(rhs, 1, IK) - 1_IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getPolyMul_RK5(lhs, rhs) result(mul)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyMul_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: lhs(:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(:)
        real(RKG)                                   :: mul(1 : min(size(lhs, 1, IK) * size(rhs, 1, IK), size(lhs, 1, IK) + size(rhs, 1, IK) - 1_IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getPolyMul_RK4(lhs, rhs) result(mul)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyMul_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: lhs(:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(:)
        real(RKG)                                   :: mul(1 : min(size(lhs, 1, IK) * size(rhs, 1, IK), size(lhs, 1, IK) + size(rhs, 1, IK) - 1_IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getPolyMul_RK3(lhs, rhs) result(mul)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyMul_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: lhs(:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(:)
        real(RKG)                                   :: mul(1 : min(size(lhs, 1, IK) * size(rhs, 1, IK), size(lhs, 1, IK) + size(rhs, 1, IK) - 1_IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getPolyMul_RK2(lhs, rhs) result(mul)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyMul_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: lhs(:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(:)
        real(RKG)                                   :: mul(1 : min(size(lhs, 1, IK) * size(rhs, 1, IK), size(lhs, 1, IK) + size(rhs, 1, IK) - 1_IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getPolyMul_RK1(lhs, rhs) result(mul)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyMul_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: lhs(:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(:)
        real(RKG)                                   :: mul(1 : min(size(lhs, 1, IK) * size(rhs, 1, IK), size(lhs, 1, IK) + size(rhs, 1, IK) - 1_IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the vector of coefficients of the polynomial resulting from the multiplication of a polynomial with another polynomial of arbitrary degrees.
    !>
    !>  \details
    !>  See the documentation of [pm_polynomial](@ref pm_polynomial) for details of the implementation.<br>
    !>
    !>  \param[out] mul     :   The output `contiguous` vector of the same type and kind as the input `lhs`, of size `size(lhs) + size(rhs) - 1`,
    !>                          containing the coefficients  (in the order of **increasing power**) of the resulting polynomial from the multiplication
    !>                          of the polynomial `lhs` of arbitrary degree with another polynomial `rhs` of arbitrary degree.<br>
    !>                          By definition, the degree of the `mul` polynomial is `size(mul) - 1`.<br>
    !>  \param[in]  lhs     :   The input `contiguous` vector of non-zero size of,<br>
    !>                          <ul>
    !>                              <li>    type `complex` of kind \CKALL, or<br>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ul>
    !>                          containing the coefficients of the left-hand-side polynomial in the multiplication, in the order of increasing power.<br>
    !>                          By definition, the degree of the `lhs` polynomial is `size(lhs) - 1`.<br>
    !>                          This means that the condition `lhs(size(lhs)) /= 0.` must hold.<br>
    !>  \param[in]  rhs     :   The input `contiguous` vector of non-zero size of the same type and kind as `lhs`,
    !>                          containing the coefficients of the right-hand-side polynomial in the multiplication, in the order of increasing power.<br>
    !>                          By definition, the degree of the `rhs` polynomial is `size(rhs) - 1`.<br>
    !>                          This means that the condition `rhs(size(rhs)) /= 0.` must hold.<br>
    !>
    !>  \interface{setPolyMul}
    !>  \code{.F90}
    !>
    !>      use pm_polynomial, only: setPolyMul
    !>
    !>      call setPolyMul(mul(1 : size(lhs) + size(rhs) - 1), lhs(:), rhs(:))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < size(rhs)` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < size(lhs)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(mul) == size(lhs) + size(rhs) - 1` must hold for the corresponding input arguments.<br>
    !>  The condition `lhs(size(lhs)) /= 0.` must hold for the corresponding input arguments.<br>
    !>  The condition `rhs(size(rhs)) /= 0.` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getPolySub](@ref pm_polynomial::getPolySub)<br>
    !>  [setPolySub](@ref pm_polynomial::setPolySub)<br>
    !>  [getPolySub](@ref pm_polynomial::getPolySub)<br>
    !>  [setPolySub](@ref pm_polynomial::setPolySub)<br>
    !>  [setPolyDiv](@ref pm_polynomial::setPolyDiv)<br>
    !>  [getPolyMul](@ref pm_polynomial::getPolyMul)<br>
    !>  [setPolyMul](@ref pm_polynomial::setPolyMul)<br>
    !>
    !>  \example{setPolyMul}
    !>  \include{lineno} example/pm_polynomial/setPolyMul/main.F90
    !>  \compilef{setPolyMul}
    !>  \output{setPolyMul}
    !>  \include{lineno} example/pm_polynomial/setPolyMul/main.out.F90
    !>
    !>  \test
    !>  [test_pm_polynomial](@ref test_pm_polynomial)
    !>
    !>  \final{setPolyMul}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 11:34 PM, August 10, 2021, Dallas, TX
    interface setPolyMul

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPolyMul_CK5(mul, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyMul_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous    :: lhs(0:)
        complex(CKG), intent(in)    , contiguous    :: rhs(0:)
        complex(CKG), intent(out)   , contiguous    :: mul(0:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPolyMul_CK4(mul, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyMul_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: lhs(0:)
        complex(CKG), intent(in)    , contiguous    :: rhs(0:)
        complex(CKG), intent(out)   , contiguous    :: mul(0:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPolyMul_CK3(mul, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyMul_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: lhs(0:)
        complex(CKG), intent(in)    , contiguous    :: rhs(0:)
        complex(CKG), intent(out)   , contiguous    :: mul(0:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPolyMul_CK2(mul, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyMul_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: lhs(0:)
        complex(CKG), intent(in)    , contiguous    :: rhs(0:)
        complex(CKG), intent(out)   , contiguous    :: mul(0:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPolyMul_CK1(mul, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyMul_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: lhs(0:)
        complex(CKG), intent(in)    , contiguous    :: rhs(0:)
        complex(CKG), intent(out)   , contiguous    :: mul(0:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPolyMul_RK5(mul, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyMul_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: lhs(0:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(0:)
        real(RKG)   , intent(out)   , contiguous    :: mul(0:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPolyMul_RK4(mul, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyMul_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: lhs(0:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(0:)
        real(RKG)   , intent(out)   , contiguous    :: mul(0:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPolyMul_RK3(mul, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyMul_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: lhs(0:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(0:)
        real(RKG)   , intent(out)   , contiguous    :: mul(0:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPolyMul_RK2(mul, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyMul_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: lhs(0:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(0:)
        real(RKG)   , intent(out)   , contiguous    :: mul(0:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPolyMul_RK1(mul, lhs, rhs)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyMul_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: lhs(0:)
        real(RKG)   , intent(in)    , contiguous    :: rhs(0:)
        real(RKG)   , intent(out)   , contiguous    :: mul(0:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the quotient and remainder of dividing a polynomial with another polynomial of arbitrary degrees.
    !>
    !>  \details
    !>  See the documentation of [pm_polynomial](@ref pm_polynomial) for details of the implementation.<br>
    !>
    !>  \param[in]  dividend    :   The input `contiguous` vector of non-zero size of,<br>
    !>                              <ul>
    !>                                  <li>    type `complex` of kind \CKALL, or<br>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ul>
    !>                              containing the coefficients of the dividend polynomial in the order of increasing power.<br>
    !>                              By definition, the degree of the dividend polynomial is `size(dividend) - 1`.<br>
    !>                              This means that the condition `dividend(size(dividend)) /= 0.` must hold.<br>
    !>  \param[in]  divisor     :   The input `contiguous` vector of non-zero size of the same type and kind as `dividend`,
    !>                              containing the coefficients of the divisor polynomial in the order of increasing power.<br>
    !>                              By definition, the degree of the divisor polynomial is `size(divisor) - 1`.<br>
    !>                              This means that the condition `divisor(size(divisor)) /= 0.` must hold.<br>
    !>  \param[out] quorem      :   The output `contiguous` vector of the same type, kind and size as the input `dividend`,
    !>                              containing the coefficients of the quotient and remainder polynomials resulting of the polynomial division.<br>
    !>                              The slice `quorem(1 : lenQuo)` contains the coefficients of the resulting quotient polynomial, in the order of increasing power.<br>
    !>                              The slice `quorem(lenQuo + 1 :)` contains the coefficients of the resulting remainder polynomial, in the order of increasing power.<br>
    !>  \param[out] lenQuo      :   The output scalar `integer` of default kind \IK containing the length of the vector of coefficients of the resulting quotient.<br>
    !>                              By definition,<br>
    !>                              <ul>
    !>                                  <li>    If the condition `dividend == divisor` holds, then `lenQuo = size(dividend)`.<br>
    !>                                  <li>    If the condition `size(dividend) < size(divisor)` holds, then `lenQuo = 0_IK`.<br>
    !>                                  <li>    If the condition `lenQuo == size(dividend)` holds, then the remainder of the division is zero,
    !>                                          as implied by the empty slice `quorem(lenQuo + 1 : size(dividend))`.<br>
    !>                              </ul>
    !>
    !>  \interface{setPolyDiv}
    !>  \code{.F90}
    !>
    !>      use pm_polynomial, only: setPolyDiv
    !>
    !>      call setPolyDiv(dividend(:), divisor(:), quorem(1:size(dividend)), lenQuo)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < size(divisor)` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < size(dividend)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(quorem) == size(dividend)` must hold for the corresponding input arguments.<br>
    !>  The condition `dividend(size(dividend)) /= 0.` must hold for the corresponding input arguments.<br>
    !>  The condition `divisor(size(divisor)) /= 0.` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getPolySub](@ref pm_polynomial::getPolySub)<br>
    !>  [setPolySub](@ref pm_polynomial::setPolySub)<br>
    !>  [getPolySub](@ref pm_polynomial::getPolySub)<br>
    !>  [setPolySub](@ref pm_polynomial::setPolySub)<br>
    !>  [setPolyDiv](@ref pm_polynomial::setPolyDiv)<br>
    !>  [getPolyMul](@ref pm_polynomial::getPolyMul)<br>
    !>  [setPolyMul](@ref pm_polynomial::setPolyMul)<br>
    !>
    !>  \example{setPolyDiv}
    !>  \include{lineno} example/pm_polynomial/setPolyDiv/main.F90
    !>  \compilef{setPolyDiv}
    !>  \output{setPolyDiv}
    !>  \include{lineno} example/pm_polynomial/setPolyDiv/main.out.F90
    !>
    !>  \test
    !>  [test_pm_polynomial](@ref test_pm_polynomial)
    !>
    !>  \final{setPolyDiv}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 11:34 PM, August 10, 2021, Dallas, TX
    interface setPolyDiv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPolyDiv_CK5(dividend, divisor, quorem, lenQuo)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiv_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous    :: dividend(:)
        complex(CKG), intent(in)    , contiguous    :: divisor(:)
        complex(CKG), intent(out)   , contiguous    :: quorem(:)
        integer(IK) , intent(out)                   :: lenQuo
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPolyDiv_CK4(dividend, divisor, quorem, lenQuo)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiv_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: dividend(:)
        complex(CKG), intent(in)    , contiguous    :: divisor(:)
        complex(CKG), intent(out)   , contiguous    :: quorem(:)
        integer(IK) , intent(out)                   :: lenQuo
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPolyDiv_CK3(dividend, divisor, quorem, lenQuo)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiv_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: dividend(:)
        complex(CKG), intent(in)    , contiguous    :: divisor(:)
        complex(CKG), intent(out)   , contiguous    :: quorem(:)
        integer(IK) , intent(out)                   :: lenQuo
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPolyDiv_CK2(dividend, divisor, quorem, lenQuo)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiv_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: dividend(:)
        complex(CKG), intent(in)    , contiguous    :: divisor(:)
        complex(CKG), intent(out)   , contiguous    :: quorem(:)
        integer(IK) , intent(out)                   :: lenQuo
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPolyDiv_CK1(dividend, divisor, quorem, lenQuo)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiv_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: dividend(:)
        complex(CKG), intent(in)    , contiguous    :: divisor(:)
        complex(CKG), intent(out)   , contiguous    :: quorem(:)
        integer(IK) , intent(out)                   :: lenQuo
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPolyDiv_RK5(dividend, divisor, quorem, lenQuo)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiv_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: dividend(:)
        real(RKG)   , intent(in)    , contiguous    :: divisor(:)
        real(RKG)   , intent(out)   , contiguous    :: quorem(:)
        integer(IK) , intent(out)                   :: lenQuo
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPolyDiv_RK4(dividend, divisor, quorem, lenQuo)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiv_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: dividend(:)
        real(RKG)   , intent(in)    , contiguous    :: divisor(:)
        real(RKG)   , intent(out)   , contiguous    :: quorem(:)
        integer(IK) , intent(out)                   :: lenQuo
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPolyDiv_RK3(dividend, divisor, quorem, lenQuo)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiv_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: dividend(:)
        real(RKG)   , intent(in)    , contiguous    :: divisor(:)
        real(RKG)   , intent(out)   , contiguous    :: quorem(:)
        integer(IK) , intent(out)                   :: lenQuo
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPolyDiv_RK2(dividend, divisor, quorem, lenQuo)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiv_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: dividend(:)
        real(RKG)   , intent(in)    , contiguous    :: divisor(:)
        real(RKG)   , intent(out)   , contiguous    :: quorem(:)
        integer(IK) , intent(out)                   :: lenQuo
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPolyDiv_RK1(dividend, divisor, quorem, lenQuo)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiv_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: dividend(:)
        real(RKG)   , intent(in)    , contiguous    :: divisor(:)
        real(RKG)   , intent(out)   , contiguous    :: quorem(:)
        integer(IK) , intent(out)                   :: lenQuo
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the vector of coefficients of the polynomial resulting from
    !>  the \f$k\f$th-order differentiation of a univariate polynomial of arbitrary degree.
    !>
    !>  \details
    !>  See the documentation of [pm_polynomial](@ref pm_polynomial) for details of the implementation.<br>
    !>
    !>  \param[in]  coef    :   The input `contiguous` vector of non-zero size of,<br>
    !>                          <ul>
    !>                              <li>    type `complex` of kind \CKALL, or<br>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ul>
    !>                          containing the coefficients (in the order of **increasing power**) of the univariate polynomial whose \f$k\f$th-order derivative must be returned.<br>
    !>                          By definition, the degree of the `coef` polynomial is `size(coef) - 1`.<br>
    !>                          This means that the condition `coef(size(coef)) /= 0.` is expected to hold (although not enforced).<br>
    !>  \param[in]  order   :   The input scalar nonnegative `integer` of default kind \IK containing the order of the derivative to compute.<br>
    !>                          (**optional**, default = `1`)
    !>
    !>  \return
    !>  `diff`              :   The output `contiguous` vector of the same type and kind as the input `coef`, of size `size(coef) - order`,
    !>                          containing the coefficients  (in the order of **increasing power**) of the resulting polynomial from the \f$k\f$th-order
    !>                          differentiation of the input polynomial with coefficients `coef` of arbitrary degree.<br>
    !>                          By definition, the degree of the `diff` polynomial is `size(diff) - order`.<br>
    !>
    !>  \interface{getPolyDiff}
    !>  \code{.F90}
    !>
    !>      use pm_polynomial, only: getPolyDiff
    !>
    !>      diff(1 : size(coef) - 1) = getPolyDiff(coef(:))
    !>      diff(1 : size(coef) - order) = getPolyDiff(coef(:), order)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= order` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getPolyDiff](@ref pm_polynomial::getPolyDiff)<br>
    !>  [setPolyDiff](@ref pm_polynomial::setPolyDiff)<br>
    !>
    !>  \example{getPolyDiff}
    !>  \include{lineno} example/pm_polynomial/getPolyDiff/main.F90
    !>  \compilef{getPolyDiff}
    !>  \output{getPolyDiff}
    !>  \include{lineno} example/pm_polynomial/getPolyDiff/main.out.F90
    !>
    !>  \test
    !>  [test_pm_polynomial](@ref test_pm_polynomial)
    !>
    !>  \final{getPolyDiff}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 11:34 PM, August 10, 2021, Dallas, TX
    interface getPolyDiff

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getPolyDiffDef_CK5(coef) result(diff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyDiffDef_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG)                                :: diff(1 : size(coef, 1, IK) - 1_IK)
    end function
#endif

#if CK4_ENABLED
    PURE module function getPolyDiffDef_CK4(coef) result(diff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyDiffDef_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG)                                :: diff(1 : size(coef, 1, IK) - 1_IK)
    end function
#endif

#if CK3_ENABLED
    PURE module function getPolyDiffDef_CK3(coef) result(diff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyDiffDef_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG)                                :: diff(1 : size(coef, 1, IK) - 1_IK)
    end function
#endif

#if CK2_ENABLED
    PURE module function getPolyDiffDef_CK2(coef) result(diff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyDiffDef_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG)                                :: diff(1 : size(coef, 1, IK) - 1_IK)
    end function
#endif

#if CK1_ENABLED
    PURE module function getPolyDiffDef_CK1(coef) result(diff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyDiffDef_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG)                                :: diff(1 : size(coef, 1, IK) - 1_IK)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getPolyDiffDef_RK5(coef) result(diff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyDiffDef_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)                                   :: diff(1 : size(coef, 1, IK) - 1_IK)
    end function
#endif

#if RK4_ENABLED
    PURE module function getPolyDiffDef_RK4(coef) result(diff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyDiffDef_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)                                   :: diff(1 : size(coef, 1, IK) - 1_IK)
    end function
#endif

#if RK3_ENABLED
    PURE module function getPolyDiffDef_RK3(coef) result(diff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyDiffDef_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)                                   :: diff(1 : size(coef, 1, IK) - 1_IK)
    end function
#endif

#if RK2_ENABLED
    PURE module function getPolyDiffDef_RK2(coef) result(diff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyDiffDef_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)                                   :: diff(1 : size(coef, 1, IK) - 1_IK)
    end function
#endif

#if RK1_ENABLED
    PURE module function getPolyDiffDef_RK1(coef) result(diff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyDiffDef_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)                                   :: diff(1 : size(coef, 1, IK) - 1_IK)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getPolyDiffOrd_CK5(coef, order) result(diff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyDiffOrd_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK) , intent(in)                    :: order
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG)                                :: diff(order : size(coef, 1, IK) - 1_IK)
    end function
#endif

#if CK4_ENABLED
    PURE module function getPolyDiffOrd_CK4(coef, order) result(diff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyDiffOrd_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK) , intent(in)                    :: order
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG)                                :: diff(order : size(coef, 1, IK) - 1_IK)
    end function
#endif

#if CK3_ENABLED
    PURE module function getPolyDiffOrd_CK3(coef, order) result(diff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyDiffOrd_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK) , intent(in)                    :: order
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG)                                :: diff(order : size(coef, 1, IK) - 1_IK)
    end function
#endif

#if CK2_ENABLED
    PURE module function getPolyDiffOrd_CK2(coef, order) result(diff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyDiffOrd_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK) , intent(in)                    :: order
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG)                                :: diff(order : size(coef, 1, IK) - 1_IK)
    end function
#endif

#if CK1_ENABLED
    PURE module function getPolyDiffOrd_CK1(coef, order) result(diff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyDiffOrd_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK) , intent(in)                    :: order
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG)                                :: diff(order : size(coef, 1, IK) - 1_IK)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getPolyDiffOrd_RK5(coef, order) result(diff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyDiffOrd_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK) , intent(in)                    :: order
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)                                   :: diff(order : size(coef, 1, IK) - 1_IK)
    end function
#endif

#if RK4_ENABLED
    PURE module function getPolyDiffOrd_RK4(coef, order) result(diff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyDiffOrd_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK) , intent(in)                    :: order
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)                                   :: diff(order : size(coef, 1, IK) - 1_IK)
    end function
#endif

#if RK3_ENABLED
    PURE module function getPolyDiffOrd_RK3(coef, order) result(diff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyDiffOrd_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK) , intent(in)                    :: order
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)                                   :: diff(order : size(coef, 1, IK) - 1_IK)
    end function
#endif

#if RK2_ENABLED
    PURE module function getPolyDiffOrd_RK2(coef, order) result(diff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyDiffOrd_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK) , intent(in)                    :: order
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)                                   :: diff(order : size(coef, 1, IK) - 1_IK)
    end function
#endif

#if RK1_ENABLED
    PURE module function getPolyDiffOrd_RK1(coef, order) result(diff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyDiffOrd_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK) , intent(in)                    :: order
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)                                   :: diff(order : size(coef, 1, IK) - 1_IK)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the vector of coefficients of the polynomial resulting from the \f$k\f$th-order differentiation of a univariate polynomial of arbitrary degree.
    !>
    !>  \details
    !>  See the documentation of [pm_polynomial](@ref pm_polynomial) for details of the implementation.<br>
    !>
    !>  \param[out] diff    :   The output `contiguous` vector of the same type and kind as the input `coef`, of size `size(coef) - order`,
    !>                          containing the coefficients  (in the order of **increasing power**) of the resulting polynomial from the \f$k\f$th-order
    !>                          differentiation of the input polynomial with coefficients `coef` of arbitrary degree.<br>
    !>                          By definition, the degree of the `diff` polynomial is `size(diff) - order`.<br>
    !>  \param[in]  coef    :   The input `contiguous` vector of non-zero size of,<br>
    !>                          <ul>
    !>                              <li>    type `complex` of kind \CKALL, or<br>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ul>
    !>                          containing the coefficients (in the order of **increasing power**) of the univariate polynomial whose \f$k\f$th-order derivative must be returned.<br>
    !>                          By definition, the degree of the `coef` polynomial is `size(coef) - 1`.<br>
    !>                          This means that the condition `coef(size(coef)) /= 0.` is expected to hold (although not enforced).<br>
    !>  \param[in]  order   :   The input scalar nonnegative `integer` of default kind \IK containing the order of the derivative to compute.<br>
    !>                          (**optional**, default = `1`)
    !>
    !>  \interface{setPolyDiff}
    !>  \code{.F90}
    !>
    !>      use pm_polynomial, only: setPolyDiff
    !>
    !>      call setPolyDiff(diff(1 : size(coef) - 1), coef(:))
    !>      call setPolyDiff(diff(1 : size(coef) - order), coef(:), order)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= order` must hold for the corresponding input arguments.<br>
    !>  The condition `size(diff) == max(1, size(coef) - order)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getPolyDiff](@ref pm_polynomial::getPolyDiff)<br>
    !>  [setPolyDiff](@ref pm_polynomial::setPolyDiff)<br>
    !>
    !>  \example{setPolyDiff}
    !>  \include{lineno} example/pm_polynomial/setPolyDiff/main.F90
    !>  \compilef{setPolyDiff}
    !>  \output{setPolyDiff}
    !>  \include{lineno} example/pm_polynomial/setPolyDiff/main.out.F90
    !>
    !>  \test
    !>  [test_pm_polynomial](@ref test_pm_polynomial)
    !>
    !>  \final{setPolyDiff}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 11:34 PM, August 10, 2021, Dallas, TX
    interface setPolyDiff

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPolyDiffDef_CK5(diff, coef)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiffDef_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG), intent(out)   , contiguous    :: diff(1:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPolyDiffDef_CK4(diff, coef)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiffDef_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG), intent(out)   , contiguous    :: diff(1:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPolyDiffDef_CK3(diff, coef)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiffDef_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG), intent(out)   , contiguous    :: diff(1:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPolyDiffDef_CK2(diff, coef)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiffDef_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG), intent(out)   , contiguous    :: diff(1:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPolyDiffDef_CK1(diff, coef)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiffDef_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG), intent(out)   , contiguous    :: diff(1:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPolyDiffDef_RK5(diff, coef)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiffDef_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)   , intent(out)   , contiguous    :: diff(1:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPolyDiffDef_RK4(diff, coef)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiffDef_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)   , intent(out)   , contiguous    :: diff(1:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPolyDiffDef_RK3(diff, coef)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiffDef_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)   , intent(out)   , contiguous    :: diff(1:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPolyDiffDef_RK2(diff, coef)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiffDef_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)   , intent(out)   , contiguous    :: diff(1:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPolyDiffDef_RK1(diff, coef)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiffDef_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)   , intent(out)   , contiguous    :: diff(1:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPolyDiffOrd_CK5(diff, coef, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiffOrd_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK) , intent(in)                    :: order
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG), intent(out)   , contiguous    :: diff(order:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPolyDiffOrd_CK4(diff, coef, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiffOrd_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK) , intent(in)                    :: order
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG), intent(out)   , contiguous    :: diff(order:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPolyDiffOrd_CK3(diff, coef, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiffOrd_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK) , intent(in)                    :: order
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG), intent(out)   , contiguous    :: diff(order:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPolyDiffOrd_CK2(diff, coef, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiffOrd_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK) , intent(in)                    :: order
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG), intent(out)   , contiguous    :: diff(order:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPolyDiffOrd_CK1(diff, coef, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiffOrd_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK) , intent(in)                    :: order
        complex(CKG), intent(in)    , contiguous    :: coef(0:)
        complex(CKG), intent(out)   , contiguous    :: diff(order:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPolyDiffOrd_RK5(diff, coef, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiffOrd_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK) , intent(in)                    :: order
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)   , intent(out)   , contiguous    :: diff(order:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPolyDiffOrd_RK4(diff, coef, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiffOrd_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK) , intent(in)                    :: order
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)   , intent(out)   , contiguous    :: diff(order:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPolyDiffOrd_RK3(diff, coef, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiffOrd_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK) , intent(in)                    :: order
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)   , intent(out)   , contiguous    :: diff(order:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPolyDiffOrd_RK2(diff, coef, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiffOrd_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK) , intent(in)                    :: order
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)   , intent(out)   , contiguous    :: diff(order:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPolyDiffOrd_RK1(diff, coef, order)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyDiffOrd_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK) , intent(in)                    :: order
        real(RKG)   , intent(in)    , contiguous    :: coef(0:)
        real(RKG)   , intent(out)   , contiguous    :: diff(order:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a string containing the polynomial expression corresponding to the input polynomial coefficients.
    !>
    !>  \param[in]  coef    :   The input `contiguous` vector of non-zero size of,<br>
    !>                          <ul>
    !>                              <li>    type `complex` of kind \CKALL, or<br>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ul>
    !>                          containing the coefficients of the polynomial in the order of increasing power.<br>
    !>
    !>  \return
    !>  `str`               :   The output `allocatable` scalar of type `character` of default kind \SK
    !>                          containing the polynomial expression corresponding to the input polynomial coefficients.
    !>
    !>  \interface{getPolyStr}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: SK
    !>      use pm_polynomial, only: getPolyStr
    !>      character(:, SK), allocatable :: str
    !>
    !>      str = getPolyStr(coef(:))
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>  The conditional impurity is caused by the call to [setResized](@ref pm_arrayResize::setResized).<br>
    !>
    !>  \see
    !>  [pm_polynomial](@ref pm_polynomial)<br>
    !>  [pm_polynomial](@ref pm_polynomial)<br>
    !>  [pm_polynomial](@ref pm_polynomial)<br>
    !>
    !>  \example{getPolyStr}
    !>  \include{lineno} example/pm_polynomial/getPolyStr/main.F90
    !>  \compilef{getPolyStr}
    !>  \output{getPolyStr}
    !>  \include{lineno} example/pm_polynomial/getPolyStr/main.out.F90
    !>
    !>  \test
    !>  [test_pm_polynomial](@ref test_pm_polynomial)
    !>
    !>  \todo
    !>  \pmed
    !>  This generic interface can be expanded to include an option for the name of the polynomial variable and operator symbols.
    !>
    !>  \todo
    !>  \plow
    !>  The runtime performance of this algorithm can be improved by calling [setStr()](@ref pm_val2str::setStr) in the implementation and adding fixed symbols (variable name, ...) progressively.<br>
    !>  In particular, the implied do-loop in the current implementation likely significantly stresses the compiler for high-degree polynomials.<br>
    !>  The significance of the improvements should be weighed in the light of the relevance of such improvements.<br>
    !>  Is this procedure going to be called frequently in numerically intensive applications? Unlikely.<br>
    !>
    !>  \final{getPolyStr}
    !>
    !>  \author
    !>  \AmirShahmoradi, Thursday 1:45 AM, August 22, 2019, Dallas, TX
    interface getPolyStr

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getPolyStrDef_CK5(coef) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyStrDef_CK5
#endif
        use pm_kind, only: SKG => SK, CKG => CK5
        complex(CKG)    , intent(in)    , contiguous    :: coef(0:)
        character(:,SKG)                , allocatable   :: str
    end function
#endif

#if CK4_ENABLED
    PURE module function getPolyStrDef_CK4(coef) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyStrDef_CK4
#endif
        use pm_kind, only: SKG => SK, CKG => CK4
        complex(CKG)    , intent(in)    , contiguous    :: coef(0:)
        character(:,SKG)                , allocatable   :: str
    end function
#endif

#if CK3_ENABLED
    PURE module function getPolyStrDef_CK3(coef) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyStrDef_CK3
#endif
        use pm_kind, only: SKG => SK, CKG => CK3
        complex(CKG)    , intent(in)    , contiguous    :: coef(0:)
        character(:,SKG)                , allocatable   :: str
    end function
#endif

#if CK2_ENABLED
    PURE module function getPolyStrDef_CK2(coef) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyStrDef_CK2
#endif
        use pm_kind, only: SKG => SK, CKG => CK2
        complex(CKG)    , intent(in)    , contiguous    :: coef(0:)
        character(:,SKG)                , allocatable   :: str
    end function
#endif

#if CK1_ENABLED
    PURE module function getPolyStrDef_CK1(coef) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyStrDef_CK1
#endif
        use pm_kind, only: SKG => SK, CKG => CK1
        complex(CKG)    , intent(in)    , contiguous    :: coef(0:)
        character(:,SKG)                , allocatable   :: str
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getPolyStrDef_RK5(coef) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyStrDef_RK5
#endif
        use pm_kind, only: SKG => SK, RKG => RK5
        real(RKG)       , intent(in)    , contiguous    :: coef(0:)
        character(:,SKG)                , allocatable   :: str
    end function
#endif

#if RK4_ENABLED
    PURE module function getPolyStrDef_RK4(coef) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyStrDef_RK4
#endif
        use pm_kind, only: SKG => SK, RKG => RK4
        real(RKG)       , intent(in)    , contiguous    :: coef(0:)
        character(:,SKG)                , allocatable   :: str
    end function
#endif

#if RK3_ENABLED
    PURE module function getPolyStrDef_RK3(coef) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyStrDef_RK3
#endif
        use pm_kind, only: SKG => SK, RKG => RK3
        real(RKG)       , intent(in)    , contiguous    :: coef(0:)
        character(:,SKG)                , allocatable   :: str
    end function
#endif

#if RK2_ENABLED
    PURE module function getPolyStrDef_RK2(coef) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyStrDef_RK2
#endif
        use pm_kind, only: SKG => SK, RKG => RK2
        real(RKG)       , intent(in)    , contiguous    :: coef(0:)
        character(:,SKG)                , allocatable   :: str
    end function
#endif

#if RK1_ENABLED
    PURE module function getPolyStrDef_RK1(coef) result(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyStrDef_RK1
#endif
        use pm_kind, only: SKG => SK, RKG => RK1
        real(RKG)       , intent(in)    , contiguous    :: coef(0:)
        character(:,SKG)                , allocatable   :: str
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is an `abstract` derived type for constructing concrete derived types to
    !>  distinguish various procedure signatures that require root-finding methods (e.g., Eigenvalue, Jenkins, Laguerre, ...).<br>
    !>
    !>  \details
    !>  This `abstract` derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users must use `parameter` objects instantiated from the concrete subclasses of this parent `abstract` derived type.<br>
    !>
    !>  \see
    !>  [sgl](@ref pm_polynomial::sgl)<br>
    !>  [eigen](@ref pm_polynomial::eigen)<br>
    !>  [jenkins](@ref pm_polynomial::jenkins)<br>
    !>  [laguerre](@ref pm_polynomial::laguerre)<br>
    !>  [sgl_type](@ref pm_polynomial::sgl_type)<br>
    !>  [eigen_type](@ref pm_polynomial::eigen_type)<br>
    !>  [jenkins_type](@ref pm_polynomial::jenkins_type)<br>
    !>  [laguerre_type](@ref pm_polynomial::laguerre_type)<br>
    !>  [method_type](@ref pm_polynomial::method_type)<br>
    !>
    !>  \final{method_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, abstract :: method_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the use of **Skowron-Gould** method of root-finding.<br>
    !>
    !>  \details
    !>  See the root-finding section in the documentation of [pm_polynomial](@ref pm_polynomial) for more details about this root-finding method.<br>
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  The approach of Skowron and Gould is outlined in their 2012 paper:
    !>  *General Complex Polynomial Root Solver and Its Further Optimization for Binary Microlenses*.
    !>
    !>  \interface{sgl_type}
    !>  \code{.F90}
    !>
    !>      use pm_polynomial, only: sgl_type
    !>      type(sgl_type), allocatable :: method
    !>
    !>      method = sgl_type(polished = polished, informed = informed)
    !>
    !>  \endcode
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [sgl](@ref pm_polynomial::sgl) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [sgl](@ref pm_polynomial::sgl)<br>
    !>  [eigen](@ref pm_polynomial::eigen)<br>
    !>  [jenkins](@ref pm_polynomial::jenkins)<br>
    !>  [laguerre](@ref pm_polynomial::laguerre)<br>
    !>  [sgl_type](@ref pm_polynomial::sgl_type)<br>
    !>  [eigen_type](@ref pm_polynomial::eigen_type)<br>
    !>  [jenkins_type](@ref pm_polynomial::jenkins_type)<br>
    !>  [laguerre_type](@ref pm_polynomial::laguerre_type)<br>
    !>  [method_type](@ref pm_polynomial::method_type)<br>
    !>
    !>  \final{sgl_type}
    !>
    !>  \author
    !>  \FatemehBagheri, August 28, 2024, 6:12 PM, NASA Goddard Space Flight Center, Washington, D.C.
    type, extends(method_type) :: sgl_type
        logical(LK) :: polished = .true._LK     !<  \public The scalar `logical` of default kind \LK. If `.true.`, the roots identified will be polished.
        logical(LK) :: reckoned = .false._LK    !<  \public The scalar `logical` of default kind \LK. If `.true.`, the search will start from the user-specified values, otherwise, from zero for all roots.
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [sgl_type](@ref pm_polynomial::sgl_type) that is exclusively used
    !>  to signify the use of **Skowron-Gould** method of root-finding within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  See the root-finding section in the documentation of [pm_polynomial](@ref pm_polynomial) for more details about this root-finding method.<br>
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [sgl](@ref pm_polynomial::sgl)<br>
    !>  [eigen](@ref pm_polynomial::eigen)<br>
    !>  [jenkins](@ref pm_polynomial::jenkins)<br>
    !>  [laguerre](@ref pm_polynomial::laguerre)<br>
    !>  [sgl_type](@ref pm_polynomial::sgl_type)<br>
    !>  [eigen_type](@ref pm_polynomial::eigen_type)<br>
    !>  [jenkins_type](@ref pm_polynomial::jenkins_type)<br>
    !>  [laguerre_type](@ref pm_polynomial::laguerre_type)<br>
    !>  [method_type](@ref pm_polynomial::method_type)<br>
    !>
    !>  \final{sgl}
    !>
    !>  \author
    !>  \FatemehBagheri, August 28, 2024, 6:12 PM, NASA Goddard Space Flight Center, Washington, D.C.
    type(sgl_type), parameter :: sgl = sgl_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: sgl
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the use of the Eigenvalue method of root-finding.<br>
    !>
    !>  \details
    !>  See the root-finding section in the documentation of [pm_polynomial](@ref pm_polynomial) for more details about this root-finding method.<br>
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [eigen](@ref pm_polynomial::eigen) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [sgl](@ref pm_polynomial::sgl)<br>
    !>  [eigen](@ref pm_polynomial::eigen)<br>
    !>  [jenkins](@ref pm_polynomial::jenkins)<br>
    !>  [laguerre](@ref pm_polynomial::laguerre)<br>
    !>  [sgl_type](@ref pm_polynomial::sgl_type)<br>
    !>  [eigen_type](@ref pm_polynomial::eigen_type)<br>
    !>  [jenkins_type](@ref pm_polynomial::jenkins_type)<br>
    !>  [laguerre_type](@ref pm_polynomial::laguerre_type)<br>
    !>  [method_type](@ref pm_polynomial::method_type)<br>
    !>
    !>  \final{eigen_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(method_type) :: eigen_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [eigen_type](@ref pm_polynomial::eigen_type) that is exclusively used
    !>  to signify the use of Eigenvalue method of root-finding within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  See the root-finding section in the documentation of [pm_polynomial](@ref pm_polynomial) for more details about this root-finding method.<br>
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [sgl](@ref pm_polynomial::sgl)<br>
    !>  [eigen](@ref pm_polynomial::eigen)<br>
    !>  [jenkins](@ref pm_polynomial::jenkins)<br>
    !>  [laguerre](@ref pm_polynomial::laguerre)<br>
    !>  [sgl_type](@ref pm_polynomial::sgl_type)<br>
    !>  [eigen_type](@ref pm_polynomial::eigen_type)<br>
    !>  [jenkins_type](@ref pm_polynomial::jenkins_type)<br>
    !>  [laguerre_type](@ref pm_polynomial::laguerre_type)<br>
    !>  [method_type](@ref pm_polynomial::method_type)<br>
    !>
    !>  \final{eigen}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(eigen_type), parameter :: eigen = eigen_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: eigen
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the use of **Jenkins-Traub** method of root-finding.<br>
    !>
    !>  \details
    !>  See the root-finding section in the documentation of [pm_polynomial](@ref pm_polynomial) for more details about this root-finding method.<br>
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [jenkins](@ref pm_polynomial::jenkins) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [sgl](@ref pm_polynomial::sgl)<br>
    !>  [eigen](@ref pm_polynomial::eigen)<br>
    !>  [jenkins](@ref pm_polynomial::jenkins)<br>
    !>  [laguerre](@ref pm_polynomial::laguerre)<br>
    !>  [sgl_type](@ref pm_polynomial::sgl_type)<br>
    !>  [eigen_type](@ref pm_polynomial::eigen_type)<br>
    !>  [jenkins_type](@ref pm_polynomial::jenkins_type)<br>
    !>  [laguerre_type](@ref pm_polynomial::laguerre_type)<br>
    !>  [method_type](@ref pm_polynomial::method_type)<br>
    !>
    !>  \final{jenkins_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(method_type) :: jenkins_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [jenkins_type](@ref pm_polynomial::jenkins_type) that is exclusively used
    !>  to signify the use of **Jenkins-Traub** method of root-finding within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  See the root-finding section in the documentation of [pm_polynomial](@ref pm_polynomial) for more details about this root-finding method.<br>
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [sgl](@ref pm_polynomial::sgl)<br>
    !>  [eigen](@ref pm_polynomial::eigen)<br>
    !>  [jenkins](@ref pm_polynomial::jenkins)<br>
    !>  [laguerre](@ref pm_polynomial::laguerre)<br>
    !>  [sgl_type](@ref pm_polynomial::sgl_type)<br>
    !>  [eigen_type](@ref pm_polynomial::eigen_type)<br>
    !>  [jenkins_type](@ref pm_polynomial::jenkins_type)<br>
    !>  [laguerre_type](@ref pm_polynomial::laguerre_type)<br>
    !>  [method_type](@ref pm_polynomial::method_type)<br>
    !>
    !>  \final{jenkins}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(jenkins_type), parameter :: jenkins = jenkins_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: jenkins
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the use of **Laguerre** method of root-finding.<br>
    !>
    !>  \details
    !>  See the root-finding section in the documentation of [pm_polynomial](@ref pm_polynomial) for more details about this root-finding method.<br>
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [laguerre](@ref pm_polynomial::laguerre) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [sgl](@ref pm_polynomial::sgl)<br>
    !>  [eigen](@ref pm_polynomial::eigen)<br>
    !>  [jenkins](@ref pm_polynomial::jenkins)<br>
    !>  [laguerre](@ref pm_polynomial::laguerre)<br>
    !>  [sgl_type](@ref pm_polynomial::sgl_type)<br>
    !>  [eigen_type](@ref pm_polynomial::eigen_type)<br>
    !>  [jenkins_type](@ref pm_polynomial::jenkins_type)<br>
    !>  [laguerre_type](@ref pm_polynomial::laguerre_type)<br>
    !>  [method_type](@ref pm_polynomial::method_type)<br>
    !>
    !>  \final{laguerre_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(method_type) :: laguerre_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [laguerre_type](@ref pm_polynomial::laguerre_type) that is exclusively used
    !>  to signify the use of **Laguerre** method of root-finding within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  See the root-finding section in the documentation of [pm_polynomial](@ref pm_polynomial) for more details about this root-finding method.<br>
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [sgl](@ref pm_polynomial::sgl)<br>
    !>  [eigen](@ref pm_polynomial::eigen)<br>
    !>  [jenkins](@ref pm_polynomial::jenkins)<br>
    !>  [laguerre](@ref pm_polynomial::laguerre)<br>
    !>  [sgl_type](@ref pm_polynomial::sgl_type)<br>
    !>  [eigen_type](@ref pm_polynomial::eigen_type)<br>
    !>  [jenkins_type](@ref pm_polynomial::jenkins_type)<br>
    !>  [laguerre_type](@ref pm_polynomial::laguerre_type)<br>
    !>  [method_type](@ref pm_polynomial::method_type)<br>
    !>
    !>  \final{laguerre}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(laguerre_type), parameter :: laguerre = laguerre_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: laguerre
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    !>  \brief
!    !>  Return the coefficients of a polynomial of arbitrary degree specified by its roots `root`.<br>
!    !>
!    !>  \details
!    !>  The strategy relies on the construction of the full polynomial by multiplying polynomial terms recursively.<br>
!    !>  Pure multiplication, however, is slow, and can be expedited by realizing that polynomial multiplication is equivalent to an FFT on the coefficients.<br>
!    !>
!    !>  \param[in]  root        :   The input `contiguous` vector of,<br>
!    !>                              <ol>
!    !>                                  <li>    type `complex` of kind \CKALL,<br>
!    !>                                  <li>    type `real` of kind \RKALL,<br>
!    !>                              </ol>
!    !>                              and the same size as the degree of the polynomial (i.e., `size(coef) - 1`),
!    !>                              containing the roots of the polynomial.<br>
!    !>  \param[out] coef        :   The output `contiguous` vector of the same type and kind as the input `root` of size `size(root) + 1`.<br>
!    !>                              On output, it contains the coefficients of the polynomial in the order of **increasing power**.<br>
!    !>                              By definition, the degree of the polynomial is `size(coef) - 1`.<br>
!    !>
!    !>  \interface{setPolyCoef}
!    !>  \code{.F90}
!    !>
!    !>      use pm_polynomial, only: setPolyCoef
!    !>
!    !>      call setPolyCoef(root(1 : degree), coef(0 : degree))
!    !>
!    !>  \endcode
!    !>
!    !>  \warning
!    !>  The condition `1 < size(coef)` must hold for the corresponding input arguments.<br>
!    !>  The condition `size(root) == size(coef) - 1` must hold for the corresponding input arguments.<br>
!    !>  \vericons
!    !>
!    !>  \impure
!    !>
!    !>  \see
!    !>  [getRoot](@ref pm_mathRoot::getRoot)<br>
!    !>  [setRoot](@ref pm_mathRoot::setRoot)<br>
!    !>  [getPolyRoot](@ref pm_polynomial::getPolyRoot)<br>
!    !>  [setPolyCoef](@ref pm_polynomial::setPolyRoot)<br>
!    !>
!    !>  \example{setPolyRoot}
!    !>  \include{lineno} example/pm_polynomial/setPolyRoot/main.F90
!    !>  \compilef{setPolyRoot}
!    !>  \output{setPolyRoot}
!    !>  \include{lineno} example/pm_polynomial/setPolyRoot/main.out.F90
!    !>
!    !>  \test
!    !>  [test_pm_polynomial](@ref test_pm_polynomial)
!    !>
!    !>  \author
!    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
!
!    ! Eigenvalue method.
!
!    interface setPolyRoot
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if CK5_ENABLED
!    module subroutine setPolyCoef_CK5_CK5(coef, root)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyCoef_CK5_CK5
!#endif
!        use pm_kind, only: CKG => CK5
!        complex(CKG)            , intent(out)   , contiguous    :: coef(:)
!        complex(CKG)            , intent(in)    , contiguous    :: root(:)
!    end subroutine
!#endif
!
!#if CK4_ENABLED
!    module subroutine setPolyCoef_CK4_CK4(coef, root)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyCoef_CK4_CK4
!#endif
!        use pm_kind, only: CKG => CK4
!        complex(CKG)            , intent(out)   , contiguous    :: coef(:)
!        complex(CKG)            , intent(in)    , contiguous    :: root(:)
!    end subroutine
!#endif
!
!#if CK3_ENABLED
!    module subroutine setPolyCoef_CK3_CK3(coef, root)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyCoef_CK3_CK3
!#endif
!        use pm_kind, only: CKG => CK3
!        complex(CKG)            , intent(out)   , contiguous    :: coef(:)
!        complex(CKG)            , intent(in)    , contiguous    :: root(:)
!    end subroutine
!#endif
!
!#if CK2_ENABLED
!    module subroutine setPolyCoef_CK2_CK2(coef, root)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyCoef_CK2_CK2
!#endif
!        use pm_kind, only: CKG => CK2
!        complex(CKG)            , intent(out)   , contiguous    :: coef(:)
!        complex(CKG)            , intent(in)    , contiguous    :: root(:)
!    end subroutine
!#endif
!
!#if CK1_ENABLED
!    module subroutine setPolyCoef_CK1_CK1(coef, root)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyCoef_CK1_CK1
!#endif
!        use pm_kind, only: CKG => CK1
!        complex(CKG)            , intent(out)   , contiguous    :: coef(:)
!        complex(CKG)            , intent(in)    , contiguous    :: root(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    module subroutine setPolyCoef_RK5_RK5(coef, root)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyCoef_RK5_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        real(RKG)               , intent(out)   , contiguous    :: coef(:)
!        real(RKG)               , intent(in)    , contiguous    :: root(:)
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    module subroutine setPolyCoef_RK4_RK4(coef, root)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyCoef_RK4_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        real(RKG)               , intent(out)   , contiguous    :: coef(:)
!        real(RKG)               , intent(in)    , contiguous    :: root(:)
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    module subroutine setPolyCoef_RK3_RK3(coef, root)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyCoef_RK3_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        real(RKG)               , intent(out)   , contiguous    :: coef(:)
!        real(RKG)               , intent(in)    , contiguous    :: root(:)
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    module subroutine setPolyCoef_RK2_RK2(coef, root)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyCoef_RK2_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        real(RKG)               , intent(out)   , contiguous    :: coef(:)
!        real(RKG)               , intent(in)    , contiguous    :: root(:)
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    module subroutine setPolyCoef_RK1_RK1(coef, root)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyCoef_RK1_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        real(RKG)               , intent(out)   , contiguous    :: coef(:)
!        real(RKG)               , intent(in)    , contiguous    :: root(:)
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the roots of a polynomial of arbitrary degree specified by its coefficients `coef`.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_polynomial](@ref pm_polynomial) for details of the root-finding methods.<br>
    !>
    !>  \param[in]  coef        :   The input `contiguous` vector of,<br>
    !>                              <ol>
    !>                                  <li>    type `complex` of kind \CKALL, or<br>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ol>
    !>                              containing the coefficients of the polynomial in the order of **increasing power**.<br>
    !>                              By definition, the degree of the polynomial is `size(coef) - 1`.<br>
    !>  \param[in]  method      :   The input scalar constant that can be any of the following:<br>
    !>                              <ol>
    !>                                  <li>    the scalar constant [eigen](@ref pm_polynomial::eigen) or a
    !>                                          constant object of type [eigen_type](@ref pm_polynomial::eigen_type)
    !>                                          implying the use of the Eigenvalue root-finding method.<br>
    !>                                  <li>    the scalar constant [jenkins](@ref pm_polynomial::jenkins) or a
    !>                                          constant object of type [jenkins_type](@ref pm_polynomial::jenkins_type)
    !>                                          implying the use of the Jenkins-Traub root-finding method.<br>
    !>                                  <li>    the scalar constant [laguerre](@ref pm_polynomial::laguerre) or a
    !>                                          constant object of type [laguerre_type](@ref pm_polynomial::laguerre_type)
    !>                                          implying the use of the Laguerre root-finding method.<br>
    !>                                  <li>    the scalar constant [sgl](@ref pm_polynomial::sgl) or a
    !>                                          constant object of type [sgl_type](@ref pm_polynomial::sgl_type)
    !>                                          implying the use of the Skowron-Gould root-finding method.<br>
    !>                                          **This method is yet to be fully implemented**.<br>
    !>                              </ol>
    !>                              **Which polynomial root-finding method should I use?**<br>
    !>                              <ol>
    !>                                  <li>    If you have a polynomial of highly varying coefficients, then the Eigenvalue method as
    !>                                          specified by [eigen_type](@ref pm_polynomial::eigen_type) is likely more reliable.<br>
    !>                                  <li>    The [Jenkins-Traub](@ref pm_polynomial) is also considered a relatively reliable
    !>                                          fast **Sure-Fire** technique for finding the roots of polynomials.<br>
    !>                              </ol>
    !>                              (**optional**, default = [eigen_type](@ref pm_polynomial::eigen_type))
    !>
    !>  \return
    !>  `root`                  :   The output `allocatable` vector of type `complex` of the same kind as the input `coef`,
    !>                              containing the roots of the polynomial specified by its coefficients `coef`.<br>
    !>                              By definition, the size of the output root is the same as the *degree* of the polynomial (i.e., `size(coef) - 1`).<br>
    !>                              However, **if the algorithm fails to find any of the roots at any stage**, it will return **only the roots found**.<br>
    !>                              An output `root` of zero size indicates the failure of the algorithm to converge or find any of the polynomial roots.<br>
    !>
    !>  \interface{getPolyRoot}
    !>  \code{.F90}
    !>
    !>      use pm_polynomial, only: getPolyRoot, eigen, jenkins, laguerre, sgl
    !>      complex(kind(coef)), allocatable :: root(:)
    !>
    !>      root = getPolyRoot(coef(:)) ! allocatable output.
    !>      root = getPolyRoot(coef(:), method) ! allocatable output.
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  All warnings and conditions associated with [setPolyRoot](@ref pm_polynomial::setPolyRoot) also apply to this generic interface.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>  The procedures of this generic interface are always `impure` when the input argument `method` is set to an object of type [eigen_type](@ref pm_polynomial::eigen_type).<br>
    !>
    !>  \see
    !>  [getRoot](@ref pm_mathRoot::getRoot)<br>
    !>  [setRoot](@ref pm_mathRoot::setRoot)<br>
    !>  [getPolyRoot](@ref pm_polynomial::getPolyRoot)<br>
    !>  [setPolyRoot](@ref pm_polynomial::setPolyRoot)<br>
    !>
    !>  \example{getPolyRoot}
    !>  \include{lineno} example/pm_polynomial/getPolyRoot/main.F90
    !>  \compilef{getPolyRoot}
    !>  \output{getPolyRoot}
    !>  \include{lineno} example/pm_polynomial/getPolyRoot/main.out.F90
    !>
    !>  \benchmarks
    !>
    !>  \benchmark{setPolyRoot, The effects of `method` on runtime efficiency}
    !>  The following program compares the runtime performance of [setPolyRoot](@ref pm_polynomial::setPolyRoot)
    !>  using different polynomial root finding algorithms.<br>
    !>
    !>  \include{lineno} benchmark/pm_polynomial/setPolyRoot/main.F90
    !>  \compilefb{setPolyRoot}
    !>  \postprocb{setPolyRoot}
    !>  \include{lineno} benchmark/pm_polynomial/setPolyRoot/main.py
    !>  \visb{setPolyRoot}
    !>  \image html benchmark/pm_polynomial/setPolyRoot/benchmark.setPolyRoot.runtime.png width=1000
    !>  \image html benchmark/pm_polynomial/setPolyRoot/benchmark.setPolyRoot.runtime.ratio.png width=1000
    !>  \image html benchmark/pm_polynomial/setPolyRoot/benchmark.setPolyRoot.root.count.png width=1000
    !>  \moralb{setPolyRoot}
    !>      -#  Among all summation algorithms, [jenkins_type](@ref pm_polynomial::jenkins_type)
    !>          appears to offer the fastest root finding algorithms.<br>
    !>      -#  The [eigen_type](@ref pm_polynomial::eigen_type) method also tends to offer excellent performance.<br>
    !>      -#  Unlike the above two, [laguerre_type](@ref pm_polynomial::laguerre_type) algorithm tends to significantly
    !>          trail behind both in performance and reliability in finding all roots of the polynomial.<br>
    !>
    !>  \test
    !>  [test_pm_polynomial](@ref test_pm_polynomial)
    !>
    !>  \final{getPolyRoot}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

    ! Default method.

    interface getPolyRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getPolyRootDef_CK5_CK5(coef) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootDef_CK5_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                            , allocatable   :: root(:)
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
    end function
#endif

#if CK4_ENABLED
    module function getPolyRootDef_CK4_CK4(coef) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootDef_CK4_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                            , allocatable   :: root(:)
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
    end function
#endif

#if CK3_ENABLED
    module function getPolyRootDef_CK3_CK3(coef) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootDef_CK3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                            , allocatable   :: root(:)
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
    end function
#endif

#if CK2_ENABLED
    module function getPolyRootDef_CK2_CK2(coef) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootDef_CK2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                            , allocatable   :: root(:)
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
    end function
#endif

#if CK1_ENABLED
    module function getPolyRootDef_CK1_CK1(coef) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootDef_CK1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                            , allocatable   :: root(:)
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getPolyRootDef_RK5_CK5(coef) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootDef_RK5_CK5
#endif
        use pm_kind, only: RKG => RK5
        complex(RKG)                            , allocatable   :: root(:)
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
    end function
#endif

#if RK4_ENABLED
    module function getPolyRootDef_RK4_CK4(coef) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootDef_RK4_CK4
#endif
        use pm_kind, only: RKG => RK4
        complex(RKG)                            , allocatable   :: root(:)
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
    end function
#endif

#if RK3_ENABLED
    module function getPolyRootDef_RK3_CK3(coef) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootDef_RK3_CK3
#endif
        use pm_kind, only: RKG => RK3
        complex(RKG)                            , allocatable   :: root(:)
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
    end function
#endif

#if RK2_ENABLED
    module function getPolyRootDef_RK2_CK2(coef) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootDef_RK2_CK2
#endif
        use pm_kind, only: RKG => RK2
        complex(RKG)                            , allocatable   :: root(:)
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
    end function
#endif

#if RK1_ENABLED
    module function getPolyRootDef_RK1_CK1(coef) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootDef_RK1_CK1
#endif
        use pm_kind, only: RKG => RK1
        complex(RKG)                            , allocatable   :: root(:)
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Eigenvalue method.

    interface getPolyRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getPolyRootEig_CK5_CK5(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootEig_CK5_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                            , allocatable   :: root(:)
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
        type(eigen_type)        , intent(in)                    :: method
    end function
#endif

#if CK4_ENABLED
    module function getPolyRootEig_CK4_CK4(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootEig_CK4_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                            , allocatable   :: root(:)
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
        type(eigen_type)        , intent(in)                    :: method
    end function
#endif

#if CK3_ENABLED
    module function getPolyRootEig_CK3_CK3(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootEig_CK3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                            , allocatable   :: root(:)
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
        type(eigen_type)        , intent(in)                    :: method
    end function
#endif

#if CK2_ENABLED
    module function getPolyRootEig_CK2_CK2(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootEig_CK2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                            , allocatable   :: root(:)
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
        type(eigen_type)        , intent(in)                    :: method
    end function
#endif

#if CK1_ENABLED
    module function getPolyRootEig_CK1_CK1(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootEig_CK1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                            , allocatable   :: root(:)
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
        type(eigen_type)        , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getPolyRootEig_RK5_CK5(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootEig_RK5_CK5
#endif
        use pm_kind, only: RKG => RK5
        complex(RKG)                            , allocatable   :: root(:)
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
        type(eigen_type)        , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    module function getPolyRootEig_RK4_CK4(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootEig_RK4_CK4
#endif
        use pm_kind, only: RKG => RK4
        complex(RKG)                            , allocatable   :: root(:)
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
        type(eigen_type)        , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    module function getPolyRootEig_RK3_CK3(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootEig_RK3_CK3
#endif
        use pm_kind, only: RKG => RK3
        complex(RKG)                            , allocatable   :: root(:)
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
        type(eigen_type)        , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    module function getPolyRootEig_RK2_CK2(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootEig_RK2_CK2
#endif
        use pm_kind, only: RKG => RK2
        complex(RKG)                            , allocatable   :: root(:)
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
        type(eigen_type)        , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    module function getPolyRootEig_RK1_CK1(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootEig_RK1_CK1
#endif
        use pm_kind, only: RKG => RK1
        complex(RKG)                            , allocatable   :: root(:)
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
        type(eigen_type)        , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Jenkins-Traub method.

    interface getPolyRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getPolyRootJen_CK5_CK5(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootJen_CK5_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                            , allocatable   :: root(:)
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
        type(jenkins_type)      , intent(in)                    :: method
    end function
#endif

#if CK4_ENABLED
    module function getPolyRootJen_CK4_CK4(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootJen_CK4_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                            , allocatable   :: root(:)
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
        type(jenkins_type)      , intent(in)                    :: method
    end function
#endif

#if CK3_ENABLED
    module function getPolyRootJen_CK3_CK3(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootJen_CK3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                            , allocatable   :: root(:)
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
        type(jenkins_type)      , intent(in)                    :: method
    end function
#endif

#if CK2_ENABLED
    module function getPolyRootJen_CK2_CK2(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootJen_CK2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                            , allocatable   :: root(:)
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
        type(jenkins_type)      , intent(in)                    :: method
    end function
#endif

#if CK1_ENABLED
    module function getPolyRootJen_CK1_CK1(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootJen_CK1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                            , allocatable   :: root(:)
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
        type(jenkins_type)      , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getPolyRootJen_RK5_CK5(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootJen_RK5_CK5
#endif
        use pm_kind, only: RKG => RK5
        complex(RKG)                            , allocatable   :: root(:)
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
        type(jenkins_type)      , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    module function getPolyRootJen_RK4_CK4(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootJen_RK4_CK4
#endif
        use pm_kind, only: RKG => RK4
        complex(RKG)                            , allocatable   :: root(:)
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
        type(jenkins_type)      , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    module function getPolyRootJen_RK3_CK3(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootJen_RK3_CK3
#endif
        use pm_kind, only: RKG => RK3
        complex(RKG)                            , allocatable   :: root(:)
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
        type(jenkins_type)      , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    module function getPolyRootJen_RK2_CK2(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootJen_RK2_CK2
#endif
        use pm_kind, only: RKG => RK2
        complex(RKG)                            , allocatable   :: root(:)
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
        type(jenkins_type)      , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    module function getPolyRootJen_RK1_CK1(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootJen_RK1_CK1
#endif
        use pm_kind, only: RKG => RK1
        complex(RKG)                            , allocatable   :: root(:)
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
        type(jenkins_type)      , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Laguerre method.

    interface getPolyRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getPolyRootLag_CK5_CK5(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootLag_CK5_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                            , allocatable   :: root(:)
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
        type(laguerre_type)     , intent(in)                    :: method
    end function
#endif

#if CK4_ENABLED
    module function getPolyRootLag_CK4_CK4(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootLag_CK4_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                            , allocatable   :: root(:)
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
        type(laguerre_type)     , intent(in)                    :: method
    end function
#endif

#if CK3_ENABLED
    module function getPolyRootLag_CK3_CK3(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootLag_CK3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                            , allocatable   :: root(:)
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
        type(laguerre_type)     , intent(in)                    :: method
    end function
#endif

#if CK2_ENABLED
    module function getPolyRootLag_CK2_CK2(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootLag_CK2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                            , allocatable   :: root(:)
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
        type(laguerre_type)     , intent(in)                    :: method
    end function
#endif

#if CK1_ENABLED
    module function getPolyRootLag_CK1_CK1(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootLag_CK1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                            , allocatable   :: root(:)
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
        type(laguerre_type)     , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getPolyRootLag_RK5_CK5(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootLag_RK5_CK5
#endif
        use pm_kind, only: RKG => RK5
        complex(RKG)                            , allocatable   :: root(:)
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
        type(laguerre_type)     , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    module function getPolyRootLag_RK4_CK4(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootLag_RK4_CK4
#endif
        use pm_kind, only: RKG => RK4
        complex(RKG)                            , allocatable   :: root(:)
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
        type(laguerre_type)     , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    module function getPolyRootLag_RK3_CK3(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootLag_RK3_CK3
#endif
        use pm_kind, only: RKG => RK3
        complex(RKG)                            , allocatable   :: root(:)
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
        type(laguerre_type)     , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    module function getPolyRootLag_RK2_CK2(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootLag_RK2_CK2
#endif
        use pm_kind, only: RKG => RK2
        complex(RKG)                            , allocatable   :: root(:)
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
        type(laguerre_type)     , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    module function getPolyRootLag_RK1_CK1(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootLag_RK1_CK1
#endif
        use pm_kind, only: RKG => RK1
        complex(RKG)                            , allocatable   :: root(:)
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
        type(laguerre_type)     , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Skowron-Gould method.

    interface getPolyRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getPolyRootSGL_CK5_CK5(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootSGL_CK5_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                            , allocatable   :: root(:)
        complex(CKG)            , intent(in)    , contiguous    :: coef(:)
        type(sgl_type)          , intent(in)                    :: method
    end function
#endif

#if CK4_ENABLED
    module function getPolyRootSGL_CK4_CK4(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootSGL_CK4_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                            , allocatable   :: root(:)
        complex(CKG)            , intent(in)    , contiguous    :: coef(:)
        type(sgl_type)          , intent(in)                    :: method
    end function
#endif

#if CK3_ENABLED
    module function getPolyRootSGL_CK3_CK3(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootSGL_CK3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                            , allocatable   :: root(:)
        complex(CKG)            , intent(in)    , contiguous    :: coef(:)
        type(sgl_type)          , intent(in)                    :: method
    end function
#endif

#if CK2_ENABLED
    module function getPolyRootSGL_CK2_CK2(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootSGL_CK2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                            , allocatable   :: root(:)
        complex(CKG)            , intent(in)    , contiguous    :: coef(:)
        type(sgl_type)          , intent(in)                    :: method
    end function
#endif

#if CK1_ENABLED
    module function getPolyRootSGL_CK1_CK1(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootSGL_CK1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                            , allocatable   :: root(:)
        complex(CKG)            , intent(in)    , contiguous    :: coef(:)
        type(sgl_type)          , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getPolyRootSGL_RK5_CK5(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootSGL_RK5_CK5
#endif
        use pm_kind, only: RKG => RK5
        complex(RKG)                            , allocatable   :: root(:)
        real(RKG)               , intent(in)    , contiguous    :: coef(:)
        type(sgl_type)          , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    module function getPolyRootSGL_RK4_CK4(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootSGL_RK4_CK4
#endif
        use pm_kind, only: RKG => RK4
        complex(RKG)                            , allocatable   :: root(:)
        real(RKG)               , intent(in)    , contiguous    :: coef(:)
        type(sgl_type)          , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    module function getPolyRootSGL_RK3_CK3(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootSGL_RK3_CK3
#endif
        use pm_kind, only: RKG => RK3
        complex(RKG)                            , allocatable   :: root(:)
        real(RKG)               , intent(in)    , contiguous    :: coef(:)
        type(sgl_type)          , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    module function getPolyRootSGL_RK2_CK2(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootSGL_RK2_CK2
#endif
        use pm_kind, only: RKG => RK2
        complex(RKG)                            , allocatable   :: root(:)
        real(RKG)               , intent(in)    , contiguous    :: coef(:)
        type(sgl_type)          , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    module function getPolyRootSGL_RK1_CK1(coef, method) result(root)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPolyRootSGL_RK1_CK1
#endif
        use pm_kind, only: RKG => RK1
        complex(RKG)                            , allocatable   :: root(:)
        real(RKG)               , intent(in)    , contiguous    :: coef(:)
        type(sgl_type)          , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !   \param[out] workspace   :   The output matrix of shape `(size(coef) - 1, size(coef) - 1)`,
    !                               of the same type and kind as `coef` that is used exclusively for scratch work when the Eigenvalue root-finding method is used.<br>
    !                               The specified storage will be used to construct the Upper Hessenberg Companion Matrix of the polynomial.<br>
    !                               Specifying this argument can lead to faster runtimes for repeated calls to this generic interface.<br>
    !                               On output, the companion matrix is completely destroyed before return.<br>
    !                               As such, its contents are useless.<br>
    !                               (**optional**, it can be present only if the input argument `method` is set to `eigen`.)

    !>  \brief
    !>  Return the roots of a polynomial of arbitrary degree specified by its coefficients `coef`.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_polynomial](@ref pm_polynomial) for details of the root-finding method.<br>
    !>
    !>  \param[inout]   root        :   The output `contiguous` vector of type `complex` of the same kind as the input
    !>                                  argument `coef` and the same size as the degree of the polynomial (i.e., `size(coef) - 1`),
    !>                                  containing the roots of the polynomial determined from the polynomial coefficients `coef`.<br>
    !>                                  If the specified input `method` argument is of type [sgl_type](@ref pm_polynomial::sgl_type),
    !>                                  the argument `root` has `intent(inout)` and the input values will be used to initialize the
    !>                                  root searching only if the condition `method%%reckoned == .true.` holds.<br>
    !>  \param[out]     count       :   The output scalar of type `integer` of default kind \IK, containing the
    !>                                  total number of roots computed and stored in the output vector slice `root(1 : count)`.<br>
    !>                                  A value of `count = size(root)` implies a successful computation of all polynomial roots.<br>
    !>                                  A value of `count = 0` implies a total failure of the algorithm in finding any roots.<br>
    !>                                  The condition `count < size(root)` occurs if either,<br>
    !>                                  <ul>
    !>                                      <li>    the algorithm fails to converge, or<br>
    !>                                      <li>    the algorithm fails to identify all roots of the polynomial, or<br>
    !>                                      <li>    the condition `coef(size(coef)) == 0.` occurs (i.e., when the coefficient of the highest power of the polynomial is zero), or<br>
    !>                                      <li>    the condition `size(coef) < 2` occurs (i.e., when the input polynomial is a constant).<br>
    !>                                  </ul>
    !>  \param[in]      coef        :   The input `contiguous` vector of,<br>
    !>                                  <ol>
    !>                                      <li>    type `complex` of kind \CKALL, or<br>
    !>                                      <li>    type `real` of kind \RKALL, <br>
    !>                                  </ol>
    !>                                  containing the coefficients of the polynomial in the order of **increasing power**.<br>
    !>                                  By definition, the degree of the polynomial is `size(coef) - 1`.<br>
    !>  \param[in]      method      :   The input scalar constant that can be any of the following:<br>
    !>                                  <ol>
    !>                                      <li>    the scalar constant [eigen](@ref pm_polynomial::eigen) or a
    !>                                              constant object of type [eigen_type](@ref pm_polynomial::eigen_type)
    !>                                              implying the use of the Eigenvalue root-finding method.<br>
    !>                                      <li>    the scalar constant [jenkins](@ref pm_polynomial::jenkins) or a
    !>                                              constant object of type [jenkins_type](@ref pm_polynomial::jenkins_type)
    !>                                              implying the use of the Jenkins-Traub root-finding method.<br>
    !>                                      <li>    the scalar constant [laguerre](@ref pm_polynomial::laguerre) or a
    !>                                              constant object of type [laguerre_type](@ref pm_polynomial::laguerre_type)
    !>                                              implying the use of the Laguerre root-finding method.<br>
    !>                                      <li>    the scalar constant [sgl](@ref pm_polynomial::sgl) or a
    !>                                              constant object of type [sgl_type](@ref pm_polynomial::sgl_type)
    !>                                              implying the use of the Skowron-Gould root-finding method.<br>
    !>                                              **This method is yet to be fully implemented**.<br>
    !>                                  </ol>
    !>                                  **Which polynomial root-finding method should I use?**<br>
    !>                                  <ol>
    !>                                      <li>    If you have a polynomial of highly varying coefficients, then the Eigenvalue method as
    !>                                              specified by [eigen_type](@ref pm_polynomial::eigen_type) is likely more reliable.<br>
    !>                                      <li>    The [Jenkins-Traub](@ref pm_polynomial) is also considered a relatively reliable
    !>                                              fast **Sure-Fire** technique for finding the roots of polynomials.<br>
    !>                                  </ol>
    !>
    !>  \interface{setPolyRoot}
    !>  \code{.F90}
    !>
    !>      use pm_polynomial, only: setPolyRoot, eigen_type, jenkins_type, laguerre_type, sgl_type
    !>
    !>      call setPolyRoot(root(1 : degree), count, coef(0 : degree), method) ! `degree` is the degree of the polynomial.
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `1 < size(coef)` must hold for the corresponding input arguments.<br>
    !>  The condition `coef(size(coef)) /= 0.` must hold for the corresponding input arguments (i.e., the coefficient of the highest-degree term must be non-zero).<br>
    !>  The condition `all(shape(workspace) == size(coef) - 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(root) == size(coef) - 1` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warning
    !>  It is crucial to keep in mind that computers use a fixed number of binary digits to represent floating-point numbers.<br>
    !>  As such polynomials of high degree can be problematic for root-finding algorithms.<br>
    !>
    !>  \impure
    !>
    !>  \remark
    !>  This generic interface combines, significantly extends, and modernizes the `SPOLZ`, `DPOLZ`, `CPOLZ`, and `ZPOLZ`
    !>  routines of the [MATH77](http://netlib.org/math) netlib public-domain library.<br>
    !>
    !>  \see
    !>  [getRoot](@ref pm_mathRoot::getRoot)<br>
    !>  [setRoot](@ref pm_mathRoot::setRoot)<br>
    !>  [getPolyRoot](@ref pm_polynomial::getPolyRoot)<br>
    !>  [setPolyRoot](@ref pm_polynomial::setPolyRoot)<br>
    !>
    !>  \example{setPolyRoot}
    !>  \include{lineno} example/pm_polynomial/setPolyRoot/main.F90
    !>  \compilef{setPolyRoot}
    !>  \output{setPolyRoot}
    !>  \include{lineno} example/pm_polynomial/setPolyRoot/main.out.F90
    !>
    !>  \benchmarks
    !>
    !>  \benchmark{setPolyRoot, The effects of `method` on runtime efficiency}
    !>  The following program compares the runtime performance of [setPolyRoot](@ref pm_polynomial::setPolyRoot)
    !>  using different polynomial root finding algorithms.<br>
    !>
    !>  \include{lineno} benchmark/pm_polynomial/setPolyRoot/main.F90
    !>  \compilefb{setPolyRoot}
    !>  \postprocb{setPolyRoot}
    !>  \include{lineno} benchmark/pm_polynomial/setPolyRoot/main.py
    !>  \visb{setPolyRoot}
    !>  \image html benchmark/pm_polynomial/setPolyRoot/benchmark.setPolyRoot.runtime.png width=1000
    !>  \image html benchmark/pm_polynomial/setPolyRoot/benchmark.setPolyRoot.runtime.ratio.png width=1000
    !>  \image html benchmark/pm_polynomial/setPolyRoot/benchmark.setPolyRoot.root.count.png width=1000
    !>  \moralb{setPolyRoot}
    !>      -#  Among all root finding algorithms, [jenkins_type](@ref pm_polynomial::jenkins_type) appears to be the fastest.<br>
    !>      -#  The [eigen_type](@ref pm_polynomial::eigen_type) method also tends to offer a comparably good performance.<br>
    !>      -#  Unlike the above two, [laguerre_type](@ref pm_polynomial::laguerre_type) algorithm tends to significantly
    !>          trail behind both in performance and reliability in finding all roots of the polynomial.<br>
    !>
    !>  \test
    !>  [test_pm_polynomial](@ref test_pm_polynomial)
    !>
    !>  \todo
    !>  \pvhigh
    !>  The current implementation relies on several local allocations,
    !>  the most important of which are called `workspace` in the current implementation.<br>
    !>  Although this generic interface used to accept a `workspace` argument, it was subsequently removed in favor of simplicity.<br>
    !>  A new set of interfaces can be added in the future to allow specification of scratch space arguments to avoid repeated local allocations.<br>
    !>  This could boost performance when this generic interface is called many times repeatedly.<br>
    !>
    !>  \todo
    !>  \pvhigh
    !>  The generic interface for the Eigenvalue method internally uses the eigenvalue computing
    !>  routine of EISPACK for upper Hessenberg matrices to compute the roots of the polynomial.<br>
    !>  This internal implementation should be substituted with the corresponding
    !>  external routine once it is implemented in the ParaMonte library.<br>
    !>
    !>  \todo
    !>  \plow
    !>  The existing implementations of the algorithms for the Jenkins-Traub method
    !>  contain relics of F77 coding style in the form of a few `goto` statements.<br>
    !>  These remaining goto statements should be carefully removed in the future.<br>
    !>
    !>  \todo
    !>  \pvhigh
    !>  The method of Skowron-Gould must be fully implemented.<br>
    !>
    !>  \todo
    !>  \pmed
    !>  A benchmark comparing the runtime performances of the `complex` vs. `real`
    !>  coefficient routines of this module should be added here.<br>
    !>
    !>  \final{setPolyRoot}
    !>  <br>
    !>  Copyright © 1996 California Institute of Technology, Pasadena, California. ALL RIGHTS RESERVED.
    !>  Based on Government Sponsored Research NAS7-03001.
    !>  Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
    !>  Redistributions of source code must retain this copyright notice, this list of conditions and the following disclaimer.
    !>  Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
    !>  disclaimer in the documentation and/or other materials provided with the distribution.
    !>  Neither the name of the California Institute of Technology (Caltech) nor the names of its contributors may be used
    !>  to endorse or promote products derived from this software without specific prior written permission.
    !>  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
    !>  BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
    !>  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    !>  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
    !>  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    !>  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    !>  For those codes indicated with a Math a la Carte copyright, the same rules apply, except without the full force of the Caltech legal team.
    !>  When citing this software we request that you also mention the names of the people who wrote the software you are using.
    !>  Designed by C. L. Lawson, JPL, May 1986. Programmed by C. L. Lawson and S. Y. Chiu, JPL, May 1986, Feb. 1987.
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    !>  \FatemehBagheri, Friday 09:51 AM, September 6, 2024, NASA Goddard Space Flight Center, Washington, D.C.

    ! Eigenvalue method.

    interface setPolyRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setPolyRootEig_CK5_CK5(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootEig_CK5_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
       !complex(CKG)            , intent(out)   , contiguous    :: workspace(:,:)
        complex(CKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(eigen_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setPolyRootEig_CK4_CK4(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootEig_CK4_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
       !complex(CKG)            , intent(out)   , contiguous    :: workspace(:,:)
        complex(CKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(eigen_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setPolyRootEig_CK3_CK3(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootEig_CK3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
       !complex(CKG)            , intent(out)   , contiguous    :: workspace(:,:)
        complex(CKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(eigen_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setPolyRootEig_CK2_CK2(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootEig_CK2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
       !complex(CKG)            , intent(out)   , contiguous    :: workspace(:,:)
        complex(CKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(eigen_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setPolyRootEig_CK1_CK1(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootEig_CK1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
       !complex(CKG)            , intent(out)   , contiguous    :: workspace(:,:)
        complex(CKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(eigen_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setPolyRootEig_RK5_CK5(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootEig_RK5_CK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
       !real(RKG)               , intent(out)   , contiguous    :: workspace(:,:)
        complex(RKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(eigen_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setPolyRootEig_RK4_CK4(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootEig_RK4_CK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
       !real(RKG)               , intent(out)   , contiguous    :: workspace(:,:)
        complex(RKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(eigen_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setPolyRootEig_RK3_CK3(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootEig_RK3_CK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
       !real(RKG)               , intent(out)   , contiguous    :: workspace(:,:)
        complex(RKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(eigen_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setPolyRootEig_RK2_CK2(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootEig_RK2_CK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
       !real(RKG)               , intent(out)   , contiguous    :: workspace(:,:)
        complex(RKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(eigen_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setPolyRootEig_RK1_CK1(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootEig_RK1_CK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
       !real(RKG)               , intent(out)   , contiguous    :: workspace(:,:)
        complex(RKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(eigen_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Jenkins-Traub method.

    interface setPolyRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setPolyRootJen_CK5_CK5(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootJen_CK5_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
       !complex(CKG)            , intent(out)   , contiguous    :: workspace(:)
        complex(CKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(jenkins_type)      , intent(in)                    :: method
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setPolyRootJen_CK4_CK4(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootJen_CK4_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
       !complex(CKG)            , intent(out)   , contiguous    :: workspace(:)
        complex(CKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(jenkins_type)      , intent(in)                    :: method
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setPolyRootJen_CK3_CK3(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootJen_CK3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
       !complex(CKG)            , intent(out)   , contiguous    :: workspace(:)
        complex(CKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(jenkins_type)      , intent(in)                    :: method
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setPolyRootJen_CK2_CK2(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootJen_CK2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
       !complex(CKG)            , intent(out)   , contiguous    :: workspace(:)
        complex(CKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(jenkins_type)      , intent(in)                    :: method
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setPolyRootJen_CK1_CK1(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootJen_CK1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
       !complex(CKG)            , intent(out)   , contiguous    :: workspace(:)
        complex(CKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(jenkins_type)      , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setPolyRootJen_RK5_CK5(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootJen_RK5_CK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
       !real(RKG)               , intent(out)   , contiguous    :: workspace(:)
        complex(RKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(jenkins_type)      , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setPolyRootJen_RK4_CK4(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootJen_RK4_CK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
       !real(RKG)               , intent(out)   , contiguous    :: workspace(:)
        complex(RKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(jenkins_type)      , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setPolyRootJen_RK3_CK3(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootJen_RK3_CK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
       !real(RKG)               , intent(out)   , contiguous    :: workspace(:)
        complex(RKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(jenkins_type)      , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setPolyRootJen_RK2_CK2(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootJen_RK2_CK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
       !real(RKG)               , intent(out)   , contiguous    :: workspace(:)
        complex(RKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(jenkins_type)      , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setPolyRootJen_RK1_CK1(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootJen_RK1_CK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
       !real(RKG)               , intent(out)   , contiguous    :: workspace(:)
        complex(RKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(jenkins_type)      , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Laguerre method.

    interface setPolyRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setPolyRootLag_CK5_CK5(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootLag_CK5_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
       !complex(CKG)            , intent(out)   , contiguous    :: workspace(:)
        complex(CKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(laguerre_type)     , intent(in)                    :: method
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setPolyRootLag_CK4_CK4(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootLag_CK4_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
       !complex(CKG)            , intent(out)   , contiguous    :: workspace(:)
        complex(CKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(laguerre_type)     , intent(in)                    :: method
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setPolyRootLag_CK3_CK3(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootLag_CK3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
       !complex(CKG)            , intent(out)   , contiguous    :: workspace(:)
        complex(CKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(laguerre_type)     , intent(in)                    :: method
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setPolyRootLag_CK2_CK2(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootLag_CK2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
       !complex(CKG)            , intent(out)   , contiguous    :: workspace(:)
        complex(CKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(laguerre_type)     , intent(in)                    :: method
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setPolyRootLag_CK1_CK1(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootLag_CK1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
       !complex(CKG)            , intent(out)   , contiguous    :: workspace(:)
        complex(CKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(laguerre_type)     , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setPolyRootLag_RK5_CK5(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootLag_RK5_CK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
       !real(RKG)               , intent(out)   , contiguous    :: workspace(:)
        complex(RKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(laguerre_type)     , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setPolyRootLag_RK4_CK4(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootLag_RK4_CK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
       !real(RKG)               , intent(out)   , contiguous    :: workspace(:)
        complex(RKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(laguerre_type)     , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setPolyRootLag_RK3_CK3(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootLag_RK3_CK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
       !real(RKG)               , intent(out)   , contiguous    :: workspace(:)
        complex(RKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(laguerre_type)     , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setPolyRootLag_RK2_CK2(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootLag_RK2_CK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
       !real(RKG)               , intent(out)   , contiguous    :: workspace(:)
        complex(RKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(laguerre_type)     , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setPolyRootLag_RK1_CK1(root, count, coef, method) !, workspace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootLag_RK1_CK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
       !real(RKG)               , intent(out)   , contiguous    :: workspace(:)
        complex(RKG)            , intent(out)   , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(laguerre_type)     , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Skowron-Gould method.

    interface setPolyRoot

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setPolyRootSGL_CK5_CK5(root, count, coef, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootSGL_CK5_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: coef(:)
        complex(CKG)            , intent(inout) , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(sgl_type)          , intent(in)                    :: method
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setPolyRootSGL_CK4_CK4(root, count, coef, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootSGL_CK4_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: coef(:)
        complex(CKG)            , intent(inout) , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(sgl_type)          , intent(in)                    :: method
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setPolyRootSGL_CK3_CK3(root, count, coef, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootSGL_CK3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: coef(:)
        complex(CKG)            , intent(inout) , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(sgl_type)          , intent(in)                    :: method
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setPolyRootSGL_CK2_CK2(root, count, coef, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootSGL_CK2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: coef(:)
        complex(CKG)            , intent(inout) , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(sgl_type)          , intent(in)                    :: method
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setPolyRootSGL_CK1_CK1(root, count, coef, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootSGL_CK1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: coef(:)
        complex(CKG)            , intent(inout) , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(sgl_type)          , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setPolyRootSGL_RK5_CK5(root, count, coef, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootSGL_RK5_CK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: coef(:)
        complex(RKG)            , intent(inout) , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(sgl_type)          , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setPolyRootSGL_RK4_CK4(root, count, coef, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootSGL_RK4_CK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: coef(:)
        complex(RKG)            , intent(inout) , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(sgl_type)          , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setPolyRootSGL_RK3_CK3(root, count, coef, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootSGL_RK3_CK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: coef(:)
        complex(RKG)            , intent(inout) , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(sgl_type)          , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setPolyRootSGL_RK2_CK2(root, count, coef, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootSGL_RK2_CK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: coef(:)
        complex(RKG)            , intent(inout) , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(sgl_type)          , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setPolyRootSGL_RK1_CK1(root, count, coef, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootSGL_RK1_CK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: coef(:)
        complex(RKG)            , intent(inout) , contiguous    :: root(:)
        integer(IK)             , intent(out)                   :: count
        type(sgl_type)          , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the polished (refined) root of a polynomial of arbitrary degree specified by its coefficients `coef`.<br>
    !>
    !>  \details
    !>  The procedures of this generic interface uses the Laguerre root
    !>  polishing method to refine the root of a polynomial of arbitrary degree.<br>
    !>  See the documentation of [pm_polynomial](@ref pm_polynomial) for details of the root-finding method.<br>
    !>
    !>  \param[inout]   root    :   The input/output scalar containing the initial guess for the polynomial root.<br>
    !>                              <ol>
    !>                                  <li>    If `coef` is of type `complex`, then `root` must be of the same type and kind as `coef`.<br>
    !>                                  <li>    If `coef` is of type `real`, then `root` can be either `real` or `complex` of the same kind as `coef`.<br>
    !>                              </ol>
    !>                              On output, this argument will contain the refined (polished) root of the polynomial \f$P\f$ such that \f$P(\ms{root}) = 0\f$.<br>
    !>  \param[in]      coef    :   The input `contiguous` vector of size at least `2`,<br>
    !>                              <ol>
    !>                                  <li>    type `complex` of kind \CKALL, or<br>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ol>
    !>                              containing the coefficients of the polynomial in the order of **increasing power**.<br>
    !>                              By definition, the degree of the polynomial is `size(coef) - 1`.<br>
    !>  \param[out]     niter   :   The output scalar of type `integer` of default kind \IK, containing
    !>                              the total number of iterations performed to polish the input `root`.<br>
    !>                              <ol>
    !>                                  <li>    A **positive non-zero** value indicates **successful convergence** within `niter` number of iterations.<br>
    !>                                  <li>    A **negative non-zero** value indicates **convergence failure** within `abs(niter)` number of iterations.<br>
    !>                                  <li>    A **zero value** indicates the **error** condition `size(coef) < 2` has occurred (i.e., when the input polynomial is a constant).<br>
    !>                              </ol>
    !>                              The value of this output argument **must always** be inspected before using the output `root`.<br>
    !>
    !>  \interface{setPolyRootPolished}
    !>  \code{.F90}
    !>
    !>      use pm_polynomial, only: setPolyRootPolished
    !>
    !>      call setPolyRootPolished(root, niter, coef(:))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `1 < size(coef)` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warning
    !>  It is crucial to keep in mind that computers use a fixed number of binary digits to represent floating-point numbers.<br>
    !>  As such polynomials of high degree can be problematic for root-finding algorithms.<br>
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getRoot](@ref pm_mathRoot::getRoot)<br>
    !>  [setRoot](@ref pm_mathRoot::setRoot)<br>
    !>  [getPolyRoot](@ref pm_polynomial::getPolyRoot)<br>
    !>  [setPolyRoot](@ref pm_polynomial::setPolyRoot)<br>
    !>
    !>  \example{setPolyRootPolished}
    !>  \include{lineno} example/pm_polynomial/setPolyRootPolished/main.F90
    !>  \compilef{setPolyRootPolished}
    !>  \output{setPolyRootPolished}
    !>  \include{lineno} example/pm_polynomial/setPolyRootPolished/main.out.F90
    !>
    !>  \test
    !>  [test_pm_polynomial](@ref test_pm_polynomial)
    !>
    !>  \todo
    !>  \pmed
    !>  An optional `reltol` may be added in the future to allow control over the error in the refined `root`.<br>
    !>
    !>  \final{setPolyRootPolished}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

    interface setPolyRootPolished

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setPolyRootPolishedLag_CK5_CK5(root, niter, coef)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootPolishedLag_CK5_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
        complex(CKG)            , intent(inout)                 :: root
        integer(IK)             , intent(out)                   :: niter
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setPolyRootPolishedLag_CK4_CK4(root, niter, coef)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootPolishedLag_CK4_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
        complex(CKG)            , intent(inout)                 :: root
        integer(IK)             , intent(out)                   :: niter
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setPolyRootPolishedLag_CK3_CK3(root, niter, coef)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootPolishedLag_CK3_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
        complex(CKG)            , intent(inout)                 :: root
        integer(IK)             , intent(out)                   :: niter
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setPolyRootPolishedLag_CK2_CK2(root, niter, coef)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootPolishedLag_CK2_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
        complex(CKG)            , intent(inout)                 :: root
        integer(IK)             , intent(out)                   :: niter
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setPolyRootPolishedLag_CK1_CK1(root, niter, coef)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootPolishedLag_CK1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: coef(0:)
        complex(CKG)            , intent(inout)                 :: root
        integer(IK)             , intent(out)                   :: niter
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPolyRootPolishedLag_RK5_CK5(root, niter, coef)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootPolishedLag_RK5_CK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
        complex(RKG)            , intent(inout)                 :: root
        integer(IK)             , intent(out)                   :: niter
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPolyRootPolishedLag_RK4_CK4(root, niter, coef)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootPolishedLag_RK4_CK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
        complex(RKG)            , intent(inout)                 :: root
        integer(IK)             , intent(out)                   :: niter
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPolyRootPolishedLag_RK3_CK3(root, niter, coef)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootPolishedLag_RK3_CK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
        complex(RKG)            , intent(inout)                 :: root
        integer(IK)             , intent(out)                   :: niter
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPolyRootPolishedLag_RK2_CK2(root, niter, coef)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootPolishedLag_RK2_CK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
        complex(RKG)            , intent(inout)                 :: root
        integer(IK)             , intent(out)                   :: niter
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPolyRootPolishedLag_RK1_CK1(root, niter, coef)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootPolishedLag_RK1_CK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
        complex(RKG)            , intent(inout)                 :: root
        integer(IK)             , intent(out)                   :: niter
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setPolyRootPolishedLag_RK5_RK5(root, niter, coef)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootPolishedLag_RK5_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
        real(RKG)               , intent(inout)                 :: root
        integer(IK)             , intent(out)                   :: niter
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setPolyRootPolishedLag_RK4_RK4(root, niter, coef)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootPolishedLag_RK4_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
        real(RKG)               , intent(inout)                 :: root
        integer(IK)             , intent(out)                   :: niter
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setPolyRootPolishedLag_RK3_RK3(root, niter, coef)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootPolishedLag_RK3_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
        real(RKG)               , intent(inout)                 :: root
        integer(IK)             , intent(out)                   :: niter
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setPolyRootPolishedLag_RK2_RK2(root, niter, coef)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootPolishedLag_RK2_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
        real(RKG)               , intent(inout)                 :: root
        integer(IK)             , intent(out)                   :: niter
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setPolyRootPolishedLag_RK1_RK1(root, niter, coef)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setPolyRootPolishedLag_RK1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: coef(0:)
        real(RKG)               , intent(inout)                 :: root
        integer(IK)             , intent(out)                   :: niter
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_polynomial ! LCOV_EXCL_LINE