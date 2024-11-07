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
!>  This module contains procedures and data types for interpolation of finite samples of data.<br>
!>
!>  \details
!>  Interpolation is a type of estimation, a method of constructing (finding) new data points based on the range of a discrete set of known data points.<br>
!>
!>  Piecewise constant interpolation
!>  ================================
!>
!>  The simplest interpolation method is to locate the nearest data value, and assign the same value.<br>
!>  In simple problems, this method is unlikely to be used, as linear interpolation (see below) is almost as easy,
!>  but in higher-dimensional multivariate interpolation, this could be a favorable choice for its speed and simplicity.<br>
!>
!>  \htmlonly
!>      <img src="pm_polation@constant.png" style="width:500px;">
!>  \endhtmlonly
!>
!>  Linear interpolation
!>  ====================
!>
!>  Linear interpolation is a method of curve fitting using linear polynomials
!>  to construct new data points within the range of a discrete set of known data points.<br>
!>
!>  If the two known points are given by the coordinates \f$(x_{0},y_{0})\f$ and \f$(x_{1},y_{1})\f$,
!>  the linear interpolant is the straight line between these points.<br>
!>  For a value \f$x\f$ in the interval \f$(x_{0},x_{1})\f$, the value \f$y\f$ along the straight line is given from the equation of slopes
!>  \f{equation}{
!>      \frac {y-y_{0}}{x-x_{0}} = \frac{y_{1}-y_{0}}{x_{1}-x_{0}} ~,
!>  \f}
!>  which can be derived geometrically from the figure on the right.<br>
!>  It is a special case of **polynomial interpolation** with \f$n = 1\f$.<br>
!>
!>  Solving this equation for \f$y\f$, which is the unknown value at \f$x\f$, gives
!>  \f{equation}{
!>      \begin{aligned} y
!>          &= y_{0} + (x - x_{0}) {\frac{y_{1} - y_{0}}{x_{1}-x_{0}}} \\
!>          &= \frac{y_{0}(x_{1}-x_{0})}{x_{1}-x_{0}} + \frac{y_{1}(x-x_{0})-y_{0}(x-x_{0})}{x_{1}-x_{0}} \\
!>          &= \frac{y_{1} x - y_{1} x_{0} - y_{0} x + y_{0} x_{0} + y_{0} x_{1} - y_{0} x_{0}}{x_{1} - x_{0}} \\
!>          &= \frac{y_{0}(x_{1} - x) + y_{1}(x - x_{0})}{x_{1} - x_{0}} ~,
!>      \end{aligned}
!>  \f}
!>  which is the formula for linear interpolation in the interval \f$(x_{0}, x_{1})\f$.<br>
!>  Outside this interval, the formula is identical to linear extrapolation.<br>
!>  This formula can also be understood as a weighted average.<br>
!>  The weights are inversely related to the distance from the end points to the unknown point; the closer point has more influence than the farther point.<br>
!>  Thus, the weights are \f$\textstyle 1 - (x - x_{0}) / (x_{1} - x_{0})\f$ and \f$\textstyle 1 - (x_{1} - x) / (x_{1} - x_{0})\f$,
!>  which are normalized distances between the unknown point and each of the end points.<br>
!>  Because these sum to \f$1\f$,
!>  \f{equation}{
!>      \begin{aligned} y
!>          &= y_{0} \left(1 - {\frac{x - x_{0}}{x_{1} - x_{0}}} \right) + y_{1} \left(1-{\frac {x_{1}-x}{x_{1}-x_{0}}}\right) \\
!>          &= y_{0} \left(1 - {\frac{x - x_{0}}{x_{1} - x_{0}}}\right) + y_{1}\left({\frac {x-x_{0}}{x_{1}-x_{0}}}\right) \\
!>          &= y_{0} \left({\frac{x_{1} - x}{x_{1} - x_{0}}}\right) + y_{1}\left({\frac{x - x_{0}}{x_{1} - x_{0}}}\right)
!>      \end{aligned}
!>  \f}
!>  yielding the formula for linear interpolation given above.<br>
!>
!>  \htmlonly
!>      <img src="pm_polation@linear.png" style="width:500px;">
!>  \endhtmlonly
!>
!>  Polynomial interpolation
!>  ========================
!>
!>  Polynomial interpolation is a generalization of linear interpolation.<br>
!>  For a sample of \f$n\f$ data points, there is exactly one polynomial of degree at most \f$n − 1\f$ going through all the data points.<br>
!>  Formally, given a set of \f$n + 1\f$ data points \f$(x_{0},y_{0}),\ldots ,(x_{n},y_{n})\f$, with no two \f$x_{j}\f$ the same,
!>  a polynomial function \f$p(x) = a_{0} + a_{1}x + \cdots + a_{n}x^{n}\f$ is said to interpolate the data if \f$p(x_{j}) = y_{j}\f$ for each \f$j\in\{0,1,\dotsc,n\}\f$.<br>
!>  There is always a unique such polynomial, commonly given explicitly as either a **Lagrange polynomial** or **Newton polynomial**.<br>
!>
!>  The polynomial interpolation error is proportional to the distance between the data points to the power \f$n\f$.<br>
!>  Furthermore, the interpolant is a polynomial and thus infinitely differentiable.<br>
!>  As such, polynomial interpolation overcomes most of the problems of linear interpolation.<br>
!>
!>  However, polynomial interpolation also has some disadvantages:<br>
!>  <ol>
!>      <li>    Calculating the interpolating polynomial is computationally expensive (see computational complexity) compared to linear interpolation.<br>
!>      <li>    Polynomial interpolation may exhibit oscillatory artifacts, especially at the end points, known as **Runge phenomenon**.<br>
!>  </ol>
!>
!>  \htmlonly
!>      <img src="pm_polation@polynomial.png" style="width:500px;">
!>  \endhtmlonly
!>
!>  The Neville algorithm
!>  ---------------------
!>
!>  The Neville algorithm, derived by the mathematician Eric Harold Neville in 1934, is used for polynomial interpolation.<br>
!>  Given \f$n + 1\f$ points \f$\{(x_i, y_i), i = 1, \ldots, n + 1\}\f$, there is a unique polynomial of degree \f$\leq n\f$ which goes through the given points.<br>
!>  The Neville algorithm approximates this unknown polynomial at a desired point \f$x\f$ **without explicitly evaluating the exact-fit polynomial coefficients**.<br>
!>  As such, the Neville algorithm may be computationally inappropriate for usage with many query points on the same input data abscissa and function values.<br>
!>  Unlike the Lagrange polynomial method, the Neville algorithm interpolation method is computationally faster and provides
!>  an estimate of the error in the output approximation.<br>
!>
!>  Spline interpolation
!>  ====================
!>
!>  Spline interpolation is a form of interpolation where the interpolant is a special type of piecewise polynomial called a **spline**.<br>
!>  That is, instead of fitting a single, high-degree polynomial to all of the values at once, spline interpolation fits low-degree polynomials
!>  to small subsets of the values, for example, fitting nine cubic polynomials between each of the pairs of ten points,
!>  instead of fitting a single degree-ten polynomial to all of them.<br>
!>  Spline interpolation is often preferred over polynomial interpolation because the
!>  interpolation error can be made small even when using low-degree polynomials for the spline.<br>
!>  Spline interpolation also avoids the problem of **Runge phenomenon**,
!>  in which oscillation can occur between points when interpolating using high-degree polynomials.<br>
!>
!>  \see
!>  [pm_polynomial](@ref pm_polynomial)<br>
!>  [bspline-fortran by Jacob Williams](https://github.com/jacobwilliams/bspline-fortran)<br>
!>  [polynomial interpolation](https://en.wikipedia.org/wiki/Polynomial_interpolation)<br>
!>
!>  \todo
!>  \phigh
!>  Routines for spline and higher-dimensional interpolation methods must be implemented.<br>
!>
!>  \test
!>  [test_pm_polation](@ref test_pm_polation)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_polation

    use pm_kind, only: SK, IK, LK, RK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_polation"

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify an interpolation using the smallest node larger than the query point within an interface of a procedure of the ParaMonte library.<br>
    !>  The name `neimean` stands for *neighbor mean*.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [neimean](@ref pm_polation::neimean)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [neimean](@ref pm_polation::neimean)<br>
    !>  [neinext](@ref pm_polation::neinext)<br>
    !>  [neiprev](@ref pm_polation::neiprev)<br>
    !>  [monopol](@ref pm_polation::monopol)<br>
    !>  [neimean_type](@ref pm_polation::neimean_type)<br>
    !>  [neinext_type](@ref pm_polation::neinext_type)<br>
    !>  [neiprev_type](@ref pm_polation::neiprev_type)<br>
    !>  [monopol_type](@ref pm_polation::monopol_type) <br>
    !>  [piwipol_type](@ref pm_polation::piwipol_type)<br>
    !>
    !>  \final{neimean_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: neimean_type; end type
    type(neimean_type), parameter :: neimean = neimean_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: neimean
#endif

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify an interpolation using the smallest node larger than the query point within an interface of a procedure of the ParaMonte library.<br>
    !>  The name `neinear` stands for *neighbor nearest*.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [neinear](@ref pm_polation::neinear)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [neimean](@ref pm_polation::neimean)<br>
    !>  [neinext](@ref pm_polation::neinext)<br>
    !>  [neiprev](@ref pm_polation::neiprev)<br>
    !>  [monopol](@ref pm_polation::monopol)<br>
    !>  [neimean_type](@ref pm_polation::neimean_type)<br>
    !>  [neinext_type](@ref pm_polation::neinext_type)<br>
    !>  [neiprev_type](@ref pm_polation::neiprev_type)<br>
    !>  [monopol_type](@ref pm_polation::monopol_type) <br>
    !>  [piwipol_type](@ref pm_polation::piwipol_type)<br>
    !>
    !>  \final{neinear_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: neinear_type; end type
    type(neinear_type), parameter :: neinear = neinear_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: neinear
#endif

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify an interpolation using the smallest node larger than the query point within an interface of a procedure of the ParaMonte library.<br>
    !>  The name `neinext` stands for *neighbor next*.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [neinext](@ref pm_polation::neinext)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [neimean](@ref pm_polation::neimean)<br>
    !>  [neinext](@ref pm_polation::neinext)<br>
    !>  [neiprev](@ref pm_polation::neiprev)<br>
    !>  [monopol](@ref pm_polation::monopol)<br>
    !>  [neimean_type](@ref pm_polation::neimean_type)<br>
    !>  [neinext_type](@ref pm_polation::neinext_type)<br>
    !>  [neiprev_type](@ref pm_polation::neiprev_type)<br>
    !>  [monopol_type](@ref pm_polation::monopol_type) <br>
    !>  [piwipol_type](@ref pm_polation::piwipol_type)<br>
    !>
    !>  \final{neinext_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: neinext_type; end type
    type(neinext_type), parameter :: neinext = neinext_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: neinext
#endif

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify an interpolation using the largest node smaller than the query point within an interface of a procedure of the ParaMonte library.<br>
    !>  The name `neiprev` stands for *neighbor previous*.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [neiprev](@ref pm_polation::neiprev)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [neimean](@ref pm_polation::neimean)<br>
    !>  [neinext](@ref pm_polation::neinext)<br>
    !>  [neiprev](@ref pm_polation::neiprev)<br>
    !>  [monopol](@ref pm_polation::monopol)<br>
    !>  [neimean_type](@ref pm_polation::neimean_type)<br>
    !>  [neinext_type](@ref pm_polation::neinext_type)<br>
    !>  [neiprev_type](@ref pm_polation::neiprev_type)<br>
    !>  [monopol_type](@ref pm_polation::monopol_type) <br>
    !>  [piwipol_type](@ref pm_polation::piwipol_type)<br>
    !>
    !>  \final{neiprev_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: neiprev_type; end type
    type(neiprev_type), parameter :: neiprev = neiprev_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: neiprev
#endif

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify an interpolation using piecewise lines within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [piwilin](@ref pm_polation::piwilin)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [neimean](@ref pm_polation::neimean)<br>
    !>  [neinext](@ref pm_polation::neinext)<br>
    !>  [neiprev](@ref pm_polation::neiprev)<br>
    !>  [monopol](@ref pm_polation::monopol)<br>
    !>  [neimean_type](@ref pm_polation::neimean_type)<br>
    !>  [neinext_type](@ref pm_polation::neinext_type)<br>
    !>  [neiprev_type](@ref pm_polation::neiprev_type)<br>
    !>  [monopol_type](@ref pm_polation::monopol_type) <br>
    !>  [piwipol_type](@ref pm_polation::piwipol_type)<br>
    !>
    !>  \final{piwilin_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: piwilin_type; end type
    type(piwilin_type), parameter :: piwilin = piwilin_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: piwilin
#endif

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify an interpolation using a single polynomial of highest degree possible given the abscissa within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [monopol](@ref pm_polation::monopol)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [neimean](@ref pm_polation::neimean)<br>
    !>  [neinext](@ref pm_polation::neinext)<br>
    !>  [neiprev](@ref pm_polation::neiprev)<br>
    !>  [monopol](@ref pm_polation::monopol)<br>
    !>  [neimean_type](@ref pm_polation::neimean_type)<br>
    !>  [neinext_type](@ref pm_polation::neinext_type)<br>
    !>  [neiprev_type](@ref pm_polation::neiprev_type)<br>
    !>  [monopol_type](@ref pm_polation::monopol_type) <br>
    !>  [piwipol_type](@ref pm_polation::piwipol_type)<br>
    !>
    !>  \final{monopol_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: monopol_type; end type
    type(monopol_type), parameter :: monopol = monopol_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: monopol
#endif

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify an interpolation using multiple **piecewise polynomial** of arbitrary degree possible given the abscissa within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \see
    !>  [neimean](@ref pm_polation::neimean)<br>
    !>  [neinext](@ref pm_polation::neinext)<br>
    !>  [neiprev](@ref pm_polation::neiprev)<br>
    !>  [monopol](@ref pm_polation::monopol)<br>
    !>  [neimean_type](@ref pm_polation::neimean_type)<br>
    !>  [neinext_type](@ref pm_polation::neinext_type)<br>
    !>  [neiprev_type](@ref pm_polation::neiprev_type)<br>
    !>  [monopol_type](@ref pm_polation::monopol_type) <br>
    !>  [piwipol_type](@ref pm_polation::piwipol_type)<br>
    !>
    !>  \final{piwipol_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: piwipol_type
        integer(IK) :: degree   !<  \public The scalar `integer` of default kind \IK, containing the degree of the piecewise polynomials.
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the approximate **polynomial interpolation** value of the input specified point `x` for the specified `method`.
    !>
    !>  \details
    !>  For polynomial interpolation, the computation relies on the **Neville algorithm**.<br>
    !>
    !>  \param[in]  method  :   The input scalar constant that can be,
    !>                          <ol>
    !>                              <li>    The scalar constant [neimean](@ref pm_polation::neimean) or a scalar object of type [neimean_type](@ref pm_polation::neimean_type)
    !>                                      implying the use of the average of the `func` values of the two nearest neighbors of the input `queryx` smaller and larger than it as the output `interp`.<br>
    !>                              <li>    The scalar constant [neinear](@ref pm_polation::neinear) or a scalar object of type [neinear_type](@ref pm_polation::neinear_type)
    !>                                      implying the use of the average of the `func` value of the nearest neighbor of the input `queryx` as the output `interp`.<br>
    !>                                      Note that the nearest neighbor in this case is measured by actual Euclidean distances of neighbors to the input `queryx`.<br>
    !>                              <li>    The scalar constant [neiprev](@ref pm_polation::neiprev) or a scalar object of type [neiprev_type](@ref pm_polation::neiprev_type)
    !>                                      implying the use of the `func` value of the largest abscissa in the input `crdx` smaller than the input `queryx` as the output `interp`.<br>
    !>                              <li>    The scalar constant [neinext](@ref pm_polation::neinext) or a scalar object of type [neinext_type](@ref pm_polation::neinext_type)
    !>                                      implying the use of the `func` value of the smallest abscissa in the input `crdx` larger than the input `queryx` as the output `interp`.<br>
    !>                              <li>    The scalar constant [piwilin](@ref pm_polation::piwilin) or a scalar object of type [piwilin_type](@ref pm_polation::piwilin_type)
    !>                                      implying the use of the **linear interpolation** of the `func` values of the two `crdx` points that bracket `queryx` as the output `interp`.<br>
    !>                                      The linear interpolation implemented in this constructor is based on the Lagrange classical formula for linear interpolation.<br>
    !>                                      Suppose an input query point \f$x\f$ falls between two nodes \f$x_i\f$ and \f$x_{i+1}\f$ with the corresponding function values
    !>                                      \f$y_i\f$ and \f$y_{i+1}\f$ and we wish to estimate the corresponding interpolated value \f$y(x)\f$, which can be computed as,
    !>                                      \f{equation*}{
    !>                                           y(x) = \frac {x - x_{i+1}} {x_i - x_{i+1}} y_i + \frac {x - x_{i}} {x_{i+1} - x_{i}} y_{i+1} ~.
    !>                                      \f}
    !>                              <li>    The scalar constant [monopol](@ref pm_polation::monopol) or a scalar object of type [monopol_type](@ref pm_polation::monopol_type)
    !>                                      implying the use of a single **polynomial interpolation** of highest degree `size(crdx) - 1` possible to all pairs of `(crdx, func)` for computing the output `interp`.<br>
    !>                                      The Neville algorithm is used to approximate the polynomial interpolation.<br>
    !>                          </ol>
    !>  \param[in]  crdx    :   The input `contiguous` vector of<br>
    !>                          <ol>
    !>                              <li>    type `real` of kind \RKALL,
    !>                          </ol>
    !>                          containing the set of abscissa in **strictly ascending order** through which the constructed polynomial must pass.<br>
    !>                          If needed, the input values for `crdx` and `sortedY` can be sorted in ascending order by calling either
    !>                          [setSorted()](@ref pm_arraySort::setSorted) or [setSorted()](@ref pm_arraySort::setSorted).
    !>  \param[in]  func    :   The input `contiguous` vector of the same type, kind, and size as `crdx`,
    !>                          containing the set of function values corresponding to the input abscissa `crdx` through which the constructed polynomial must pass.
    !>  \param[in]  queryx  :   The input scalar or vector of the same type and kind as `crdx`, containing the abscissa (x-value) of the queryx point.
    !>
    !>  \return
    !>
    !>  `interp`            :   The output object of the same type, kind, rank, and shape as `queryx`,
    !>                          containing the (approximate) interpolated y-value(s) corresponding to the `queryx` point(s).
    !>
    !>  \interface{getInterp}
    !>  \code{.F90}
    !>
    !>      use pm_polation, only: getInterp
    !>
    !>      interp = getInterp(method, crdx(1:nsam), func(1:nsam), queryx)
    !>      interp(1:nque) = getInterp(method, crdx(1:nsam), func(1:nsam), queryx(1:nque))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(crdx) == size(func)` must hold for the corresponding input arguments.<br>
    !>  The condition `all([minval(crdx, 1) <= queryx .and. queryx <= maxval(crdx, 1))` must hold for the corresponding input arguments.<br>
    !>  The condition `isAscendingAll(crdx) .or. same_type_as(method, monopol_type())` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getExtrap](@ref pm_polation::getExtrap)<br>
    !>  [setExtrap](@ref pm_polation::setExtrap)<br>
    !>  [getInterp](@ref pm_polation::getInterp)<br>
    !>  [setInterp](@ref pm_polation::setInterp)<br>
    !>  [pm_sampleQuan](@ref pm_sampleQuan)<br>
    !>  [pm_arraySort](@ref pm_arraySort)<br>
    !>  [pm_quadRomb](@ref pm_quadRomb)<br>
    !>
    !>  \example{getInterp}
    !>  \include{lineno} example/pm_polation/getInterp/main.F90
    !>  \compilef{getInterp}
    !>  \output{getInterp}
    !>  \include{lineno} example/pm_polation/getInterp/main.out.F90
    !>  \postproc{getInterp}
    !>  \include{lineno} example/pm_polation/getInterp/main.py
    !>  \vis{getInterp}
    !>  \image html pm_polation/getInterp/getInterp.neimean.interp.png width=700
    !>  \image html pm_polation/getInterp/getInterp.neinear.interp.png width=700
    !>  \image html pm_polation/getInterp/getInterp.neinext.interp.png width=700
    !>  \image html pm_polation/getInterp/getInterp.neiprev.interp.png width=700
    !>  \image html pm_polation/getInterp/getInterp.piwilin.interp.png width=700
    !>  \image html pm_polation/getInterp/getInterp.monopol.interp.png width=700
    !>  \image html pm_polation/getInterp/getInterp.rungeEffect.interp.png width=700
    !>
    !>  \test
    !>  [test_pm_polation](@ref test_pm_polation)
    !>
    !>  \final{getInterp}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! monopol MNPLD

    interface getInterp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getInterpMNPLD_ND1_QD0_RK5(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpMNPLD_ND1_QD0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(monopol_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getInterpMNPLD_ND1_QD0_RK4(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpMNPLD_ND1_QD0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(monopol_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getInterpMNPLD_ND1_QD0_RK3(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpMNPLD_ND1_QD0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(monopol_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getInterpMNPLD_ND1_QD0_RK2(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpMNPLD_ND1_QD0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(monopol_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getInterpMNPLD_ND1_QD0_RK1(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpMNPLD_ND1_QD0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(monopol_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getInterpMNPLD_ND1_QD1_RK5(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpMNPLD_ND1_QD1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(monopol_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getInterpMNPLD_ND1_QD1_RK4(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpMNPLD_ND1_QD1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(monopol_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getInterpMNPLD_ND1_QD1_RK3(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpMNPLD_ND1_QD1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(monopol_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getInterpMNPLD_ND1_QD1_RK2(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpMNPLD_ND1_QD1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(monopol_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getInterpMNPLD_ND1_QD1_RK1(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpMNPLD_ND1_QD1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(monopol_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! piwilin

    interface getInterp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getInterpPWLN_ND1_QD0_RK5(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpPWLN_ND1_QD0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getInterpPWLN_ND1_QD0_RK4(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpPWLN_ND1_QD0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getInterpPWLN_ND1_QD0_RK3(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpPWLN_ND1_QD0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getInterpPWLN_ND1_QD0_RK2(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpPWLN_ND1_QD0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getInterpPWLN_ND1_QD0_RK1(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpPWLN_ND1_QD0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getInterpPWLN_ND1_QD1_RK5(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpPWLN_ND1_QD1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getInterpPWLN_ND1_QD1_RK4(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpPWLN_ND1_QD1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getInterpPWLN_ND1_QD1_RK3(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpPWLN_ND1_QD1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getInterpPWLN_ND1_QD1_RK2(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpPWLN_ND1_QD1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getInterpPWLN_ND1_QD1_RK1(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpPWLN_ND1_QD1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! neimean

    interface getInterp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getInterpMEAN_ND1_QD0_RK5(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpMEAN_ND1_QD0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getInterpMEAN_ND1_QD0_RK4(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpMEAN_ND1_QD0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getInterpMEAN_ND1_QD0_RK3(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpMEAN_ND1_QD0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getInterpMEAN_ND1_QD0_RK2(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpMEAN_ND1_QD0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getInterpMEAN_ND1_QD0_RK1(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpMEAN_ND1_QD0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getInterpMEAN_ND1_QD1_RK5(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpMEAN_ND1_QD1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getInterpMEAN_ND1_QD1_RK4(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpMEAN_ND1_QD1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getInterpMEAN_ND1_QD1_RK3(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpMEAN_ND1_QD1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getInterpMEAN_ND1_QD1_RK2(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpMEAN_ND1_QD1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getInterpMEAN_ND1_QD1_RK1(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpMEAN_ND1_QD1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! neinear

    interface getInterp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getInterpNEAR_ND1_QD0_RK5(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpNEAR_ND1_QD0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getInterpNEAR_ND1_QD0_RK4(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpNEAR_ND1_QD0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getInterpNEAR_ND1_QD0_RK3(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpNEAR_ND1_QD0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getInterpNEAR_ND1_QD0_RK2(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpNEAR_ND1_QD0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getInterpNEAR_ND1_QD0_RK1(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpNEAR_ND1_QD0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getInterpNEAR_ND1_QD1_RK5(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpNEAR_ND1_QD1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getInterpNEAR_ND1_QD1_RK4(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpNEAR_ND1_QD1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getInterpNEAR_ND1_QD1_RK3(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpNEAR_ND1_QD1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getInterpNEAR_ND1_QD1_RK2(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpNEAR_ND1_QD1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getInterpNEAR_ND1_QD1_RK1(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpNEAR_ND1_QD1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! neinext

    interface getInterp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getInterpNEXT_ND1_QD0_RK5(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpNEXT_ND1_QD0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getInterpNEXT_ND1_QD0_RK4(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpNEXT_ND1_QD0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getInterpNEXT_ND1_QD0_RK3(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpNEXT_ND1_QD0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getInterpNEXT_ND1_QD0_RK2(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpNEXT_ND1_QD0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getInterpNEXT_ND1_QD0_RK1(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpNEXT_ND1_QD0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getInterpNEXT_ND1_QD1_RK5(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpNEXT_ND1_QD1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getInterpNEXT_ND1_QD1_RK4(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpNEXT_ND1_QD1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getInterpNEXT_ND1_QD1_RK3(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpNEXT_ND1_QD1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getInterpNEXT_ND1_QD1_RK2(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpNEXT_ND1_QD1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getInterpNEXT_ND1_QD1_RK1(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpNEXT_ND1_QD1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! neiprev

    interface getInterp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getInterpPREV_ND1_QD0_RK5(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpPREV_ND1_QD0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getInterpPREV_ND1_QD0_RK4(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpPREV_ND1_QD0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getInterpPREV_ND1_QD0_RK3(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpPREV_ND1_QD0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getInterpPREV_ND1_QD0_RK2(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpPREV_ND1_QD0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getInterpPREV_ND1_QD0_RK1(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpPREV_ND1_QD0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: interp
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getInterpPREV_ND1_QD1_RK5(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpPREV_ND1_QD1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getInterpPREV_ND1_QD1_RK4(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpPREV_ND1_QD1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getInterpPREV_ND1_QD1_RK3(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpPREV_ND1_QD1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getInterpPREV_ND1_QD1_RK2(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpPREV_ND1_QD1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getInterpPREV_ND1_QD1_RK1(method, crdx, func, queryx) result(interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInterpPREV_ND1_QD1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: interp(size(queryx, 1, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the approximate **polynomial interpolation** value of the input specified point `x` for the specified `method`.
    !>
    !>  \details
    !>  For polynomial interpolation, the computation relies on the **Neville algorithm**.<br>
    !>
    !>  \param[in]  method  :   The input scalar constant that can be,
    !>                          <ol>
    !>                              <li>    The scalar constant [neimean](@ref pm_polation::neimean) or a scalar object of type [neimean_type](@ref pm_polation::neimean_type)
    !>                                      implying the use of the average of the `func` values of the two nearest neighbors of the input `queryx` smaller and larger than it as the output `interp`.<br>
    !>                              <li>    The scalar constant [neinear](@ref pm_polation::neinear) or a scalar object of type [neinear_type](@ref pm_polation::neinear_type)
    !>                                      implying the use of the average of the `func` value of the `neinear` nearest neighbor of the input `queryx` as the output `interp`.<br>
    !>                                      Note that the nearest neighbor in this case is measured by actual Euclidean distances of neighbors to the input `queryx`.<br>
    !>                              <li>    The scalar constant [neiprev](@ref pm_polation::neiprev) or a scalar object of type [neiprev_type](@ref pm_polation::neiprev_type)
    !>                                      implying the use of the `func` value of the largest abscissa in the input `crdx` smaller than the input `queryx` as the output `interp`.<br>
    !>                              <li>    The scalar constant [neinext](@ref pm_polation::neinext) or a scalar object of type [neinext_type](@ref pm_polation::neinext_type)
    !>                                      implying the use of the `func` value of the smallest abscissa in the input `crdx` larger than the input `queryx` as the output `interp`.<br>
    !>                              <li>    The scalar constant [piwilin](@ref pm_polation::piwilin) or a scalar object of type [piwilin_type](@ref pm_polation::piwilin_type)
    !>                                      implying the use of the **linear interpolation** of the `func` values of the two `crdx` points that bracket `queryx` as the output `interp`.<br>
    !>                                      The linear interpolation implemented in this constructor is based on the Lagrange classical formula for linear interpolation.<br>
    !>                                      Suppose an input query point \f$x\f$ falls between two nodes \f$x_i\f$ and \f$x_{i+1}\f$ with the corresponding function values
    !>                                      \f$y_i\f$ and \f$y_{i+1}\f$ and we wish to estimate the corresponding interpolated value \f$y(x)\f$, which can be computed as,
    !>                                      \f{equation*}{
    !>                                           y(x) = \frac {x - x_{i+1}} {x_i - x_{i+1}} y_i + \frac {x - x_{i}} {x_{i+1} - x_{i}} y_{i+1} ~.
    !>                                      \f}
    !>                              <li>    The scalar constant [monopol](@ref pm_polation::monopol) or a scalar object of type [monopol_type](@ref pm_polation::monopol_type)
    !>                                      implying the use of a single **polynomial interpolation** of highest degree `size(crdx) - 1` possible to all pairs of `(crdx, func)` for computing the output `interp`.<br>
    !>                                      The Neville algorithm is used to approximate the polynomial interpolation.<br>
    !>                          </ol>
    !>  \param[in]  crdx    :   The input `contiguous` vector of<br>
    !>                          <ol>
    !>                              <li>    type `real` of kind \RKALL,
    !>                          </ol>
    !>                          containing the set of abscissa in **strictly ascending order** through which the constructed polynomial must pass.<br>
    !>                          If needed, the input values for `crdx` and `sortedY` can be sorted in ascending order by calling either
    !>                          [setSorted()](@ref pm_arraySort::setSorted) or [setSorted()](@ref pm_arraySort::setSorted).
    !>  \param[in]  func    :   The input `contiguous` vector of the same type, kind, and size as `crdx`,
    !>                          containing the set of function values corresponding to the input abscissa `crdx` through which the constructed polynomial must pass.
    !>  \param[in]  queryx  :   The input scalar or vector of the same type and kind as `crdx`, containing the abscissa (x-value) of the queryx point.
    !>  \param[out] interp  :   The output object of the same type, kind, rank, and shape as `queryx`,
    !>                          containing the (approximate) interpolated y-value(s) corresponding to the `queryx` point(s).
    !>  \param[out] relerr  :   The output scalar of the same type and kind as `crdx`, containing the estimated error in the output `interp`.
    !>                          (**optional**. It can be present **if and only** the input argument `method` is set to
    !>                          [monopol](@ref pm_polation::monopol) or a scalar object of type [monopol_type](@ref pm_polation::monopol_type) .)
    !>
    !>  \interface{setInterp}
    !>  \code{.F90}
    !>
    !>      use pm_polation, only: setInterp
    !>
    !>      call setInterp(method, crdx(1:nsam), func(1:nsam), queryx, interp)
    !>      call setInterp(method, crdx(1:nsam), func(1:nsam), queryx(1:nque), interp(1:nque))
    !>
    !>      call setInterp(method, crdx(1:nsam), func(1:nsam), queryx, interp, relerr) ! method = monopol
    !>      call setInterp(method, crdx(1:nsam), func(1:nsam), queryx(1:nque), interp(1:nque), relerr(1:nque)) ! method = monopol
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(crdx) == size(func)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(shape(queryx) == shape(interp))` must hold for the corresponding input arguments.<br>
    !>  The condition `all(shape(queryx) == shape(relerr))` must hold for the corresponding input arguments.<br>
    !>  The condition `all([minval(crdx, 1) <= queryx .and. queryx <= maxval(crdx, 1))` must hold for the corresponding input arguments.<br>
    !>  The condition `isAscendingAll(crdx) .or. same_type_as(method, monopol_type())` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getExtrap](@ref pm_polation::getExtrap)<br>
    !>  [setExtrap](@ref pm_polation::setExtrap)<br>
    !>  [getInterp](@ref pm_polation::getInterp)<br>
    !>  [setInterp](@ref pm_polation::setInterp)<br>
    !>  [pm_sampleQuan](@ref pm_sampleQuan)<br>
    !>  [pm_arraySort](@ref pm_arraySort)<br>
    !>  [pm_quadRomb](@ref pm_quadRomb)<br>
    !>
    !>  \example{setInterp}
    !>  \include{lineno} example/pm_polation/setInterp/main.F90
    !>  \compilef{setInterp}
    !>  \output{setInterp}
    !>  \include{lineno} example/pm_polation/setInterp/main.out.F90
    !>  \postproc{setInterp}
    !>  \include{lineno} example/pm_polation/setInterp/main.py
    !>  \vis{setInterp}
    !>  \image html pm_polation/setInterp/setInterp.neimean.interp.png width=700
    !>  \image html pm_polation/setInterp/setInterp.neinear.interp.png width=700
    !>  \image html pm_polation/setInterp/setInterp.neinext.interp.png width=700
    !>  \image html pm_polation/setInterp/setInterp.neiprev.interp.png width=700
    !>  \image html pm_polation/setInterp/setInterp.piwilin.interp.png width=700
    !>  \image html pm_polation/setInterp/setInterp.monopol.interp.png width=700
    !>  \image html pm_polation/setInterp/setInterp.rungeEffect.interp.png width=700
    !>
    !>  \test
    !>  [test_pm_polation](@ref test_pm_polation)
    !>
    !>  \final{setInterp}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! monopol MNPLD

    interface setInterp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setInterpMNPLD_ND1_QD0_RK5(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMNPLD_ND1_QD0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setInterpMNPLD_ND1_QD0_RK4(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMNPLD_ND1_QD0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setInterpMNPLD_ND1_QD0_RK3(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMNPLD_ND1_QD0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setInterpMNPLD_ND1_QD0_RK2(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMNPLD_ND1_QD0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setInterpMNPLD_ND1_QD0_RK1(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMNPLD_ND1_QD0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setInterpMNPLD_ND1_QD1_RK5(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMNPLD_ND1_QD1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setInterpMNPLD_ND1_QD1_RK4(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMNPLD_ND1_QD1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setInterpMNPLD_ND1_QD1_RK3(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMNPLD_ND1_QD1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setInterpMNPLD_ND1_QD1_RK2(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMNPLD_ND1_QD1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setInterpMNPLD_ND1_QD1_RK1(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMNPLD_ND1_QD1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! monopol MNPLE

    interface setInterp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setInterpMNPLE_ND1_QD0_RK5(method, crdx, func, queryx, interp, relerr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMNPLE_ND1_QD0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp, relerr
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setInterpMNPLE_ND1_QD0_RK4(method, crdx, func, queryx, interp, relerr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMNPLE_ND1_QD0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp, relerr
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setInterpMNPLE_ND1_QD0_RK3(method, crdx, func, queryx, interp, relerr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMNPLE_ND1_QD0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp, relerr
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setInterpMNPLE_ND1_QD0_RK2(method, crdx, func, queryx, interp, relerr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMNPLE_ND1_QD0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp, relerr
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setInterpMNPLE_ND1_QD0_RK1(method, crdx, func, queryx, interp, relerr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMNPLE_ND1_QD0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp, relerr
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setInterpMNPLE_ND1_QD1_RK5(method, crdx, func, queryx, interp, relerr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMNPLE_ND1_QD1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:), relerr(:)
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setInterpMNPLE_ND1_QD1_RK4(method, crdx, func, queryx, interp, relerr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMNPLE_ND1_QD1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:), relerr(:)
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setInterpMNPLE_ND1_QD1_RK3(method, crdx, func, queryx, interp, relerr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMNPLE_ND1_QD1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:), relerr(:)
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setInterpMNPLE_ND1_QD1_RK2(method, crdx, func, queryx, interp, relerr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMNPLE_ND1_QD1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:), relerr(:)
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setInterpMNPLE_ND1_QD1_RK1(method, crdx, func, queryx, interp, relerr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMNPLE_ND1_QD1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:), relerr(:)
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! piwilin

    interface setInterp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setInterpPWLN_ND1_QD0_RK5(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpPWLN_ND1_QD0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(piwilin_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setInterpPWLN_ND1_QD0_RK4(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpPWLN_ND1_QD0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(piwilin_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setInterpPWLN_ND1_QD0_RK3(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpPWLN_ND1_QD0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(piwilin_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setInterpPWLN_ND1_QD0_RK2(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpPWLN_ND1_QD0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(piwilin_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setInterpPWLN_ND1_QD0_RK1(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpPWLN_ND1_QD0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(piwilin_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setInterpPWLN_ND1_QD1_RK5(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpPWLN_ND1_QD1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(piwilin_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setInterpPWLN_ND1_QD1_RK4(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpPWLN_ND1_QD1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(piwilin_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setInterpPWLN_ND1_QD1_RK3(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpPWLN_ND1_QD1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(piwilin_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setInterpPWLN_ND1_QD1_RK2(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpPWLN_ND1_QD1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(piwilin_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setInterpPWLN_ND1_QD1_RK1(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpPWLN_ND1_QD1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(piwilin_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! neimean

    interface setInterp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setInterpMEAN_ND1_QD0_RK5(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMEAN_ND1_QD0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(neimean_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setInterpMEAN_ND1_QD0_RK4(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMEAN_ND1_QD0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(neimean_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setInterpMEAN_ND1_QD0_RK3(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMEAN_ND1_QD0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(neimean_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setInterpMEAN_ND1_QD0_RK2(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMEAN_ND1_QD0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(neimean_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setInterpMEAN_ND1_QD0_RK1(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMEAN_ND1_QD0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(neimean_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setInterpMEAN_ND1_QD1_RK5(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMEAN_ND1_QD1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(neimean_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setInterpMEAN_ND1_QD1_RK4(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMEAN_ND1_QD1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(neimean_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setInterpMEAN_ND1_QD1_RK3(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMEAN_ND1_QD1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(neimean_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setInterpMEAN_ND1_QD1_RK2(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMEAN_ND1_QD1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(neimean_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setInterpMEAN_ND1_QD1_RK1(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpMEAN_ND1_QD1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(neimean_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! neinear

    interface setInterp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setInterpNEAR_ND1_QD0_RK5(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpNEAR_ND1_QD0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(neinear_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setInterpNEAR_ND1_QD0_RK4(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpNEAR_ND1_QD0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(neinear_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setInterpNEAR_ND1_QD0_RK3(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpNEAR_ND1_QD0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(neinear_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setInterpNEAR_ND1_QD0_RK2(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpNEAR_ND1_QD0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(neinear_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setInterpNEAR_ND1_QD0_RK1(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpNEAR_ND1_QD0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(neinear_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setInterpNEAR_ND1_QD1_RK5(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpNEAR_ND1_QD1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(neinear_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setInterpNEAR_ND1_QD1_RK4(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpNEAR_ND1_QD1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(neinear_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setInterpNEAR_ND1_QD1_RK3(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpNEAR_ND1_QD1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(neinear_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setInterpNEAR_ND1_QD1_RK2(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpNEAR_ND1_QD1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(neinear_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setInterpNEAR_ND1_QD1_RK1(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpNEAR_ND1_QD1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(neinear_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! neinext

    interface setInterp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setInterpNEXT_ND1_QD0_RK5(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpNEXT_ND1_QD0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(neinext_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setInterpNEXT_ND1_QD0_RK4(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpNEXT_ND1_QD0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(neinext_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setInterpNEXT_ND1_QD0_RK3(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpNEXT_ND1_QD0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(neinext_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setInterpNEXT_ND1_QD0_RK2(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpNEXT_ND1_QD0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(neinext_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setInterpNEXT_ND1_QD0_RK1(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpNEXT_ND1_QD0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(neinext_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setInterpNEXT_ND1_QD1_RK5(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpNEXT_ND1_QD1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(neinext_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setInterpNEXT_ND1_QD1_RK4(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpNEXT_ND1_QD1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(neinext_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setInterpNEXT_ND1_QD1_RK3(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpNEXT_ND1_QD1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(neinext_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setInterpNEXT_ND1_QD1_RK2(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpNEXT_ND1_QD1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(neinext_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setInterpNEXT_ND1_QD1_RK1(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpNEXT_ND1_QD1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(neinext_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! neiprev

    interface setInterp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setInterpPREV_ND1_QD0_RK5(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpPREV_ND1_QD0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(neiprev_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setInterpPREV_ND1_QD0_RK4(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpPREV_ND1_QD0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(neiprev_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setInterpPREV_ND1_QD0_RK3(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpPREV_ND1_QD0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(neiprev_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setInterpPREV_ND1_QD0_RK2(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpPREV_ND1_QD0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(neiprev_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setInterpPREV_ND1_QD0_RK1(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpPREV_ND1_QD0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: interp
        type(neiprev_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setInterpPREV_ND1_QD1_RK5(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpPREV_ND1_QD1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(neiprev_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setInterpPREV_ND1_QD1_RK4(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpPREV_ND1_QD1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(neiprev_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setInterpPREV_ND1_QD1_RK3(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpPREV_ND1_QD1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(neiprev_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setInterpPREV_ND1_QD1_RK2(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpPREV_ND1_QD1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(neiprev_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setInterpPREV_ND1_QD1_RK1(method, crdx, func, queryx, interp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterpPREV_ND1_QD1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: interp(:)
        type(neiprev_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the approximate **polynomial interpolation/extrapolation** value of the input specified point `x` for the specified `method`.
    !>
    !>  \details
    !>  For polynomial interpolation/extrapolation, the computation relies on the **Neville algorithm**.<br>
    !>  The extrapolation is done as if the out-of-bound query point is within the boundary (first or last)
    !>  interpolation segment specified by the input `(crdx, func)` pairs of values.<br>
    !>
    !>  \param[in]  method  :   The input scalar constant that can be,
    !>                          <ol>
    !>                              <li>    The scalar constant [neimean](@ref pm_polation::neimean) or a scalar object of type [neimean_type](@ref pm_polation::neimean_type)
    !>                                      implying the use of the average of the `func` values of the two nearest neighbors of the input `queryx` smaller and larger than it as the output `extrap`.<br>
    !>                              <li>    The scalar constant [neinear](@ref pm_polation::neinear) or a scalar object of type [neinear_type](@ref pm_polation::neinear_type)
    !>                                      implying the use of the average of the `func` value of the `neinear` nearest neighbor of the input `queryx` as the output `extrap`.<br>
    !>                                      Note that the nearest neighbor in this case is measured by actual Euclidean distances of neighbors to the input `queryx`.<br>
    !>                              <li>    The scalar constant [neiprev](@ref pm_polation::neiprev) or a scalar object of type [neiprev_type](@ref pm_polation::neiprev_type)
    !>                                      implying the use of the `func` value of the largest abscissa in the input `crdx` smaller than the input `queryx` as the output `extrap`.<br>
    !>                              <li>    The scalar constant [neinext](@ref pm_polation::neinext) or a scalar object of type [neinext_type](@ref pm_polation::neinext_type)
    !>                                      implying the use of the `func` value of the smallest abscissa in the input `crdx` larger than the input `queryx` as the output `extrap`.<br>
    !>                              <li>    The scalar constant [piwilin](@ref pm_polation::piwilin) or a scalar object of type [piwilin_type](@ref pm_polation::piwilin_type)
    !>                                      implying the use of the **linear interpolation/extrapolation** of the `func` values of the two `crdx` points that bracket `queryx` as the output `extrap`.<br>
    !>                                      The linear interpolation/extrapolation implemented in this constructor is based on the Lagrange classical formula for linear interpolation/extrapolation.<br>
    !>                                      Suppose an input query point \f$x\f$ falls between two nodes \f$x_i\f$ and \f$x_{i+1}\f$ with the corresponding function values
    !>                                      \f$y_i\f$ and \f$y_{i+1}\f$ and we wish to estimate the corresponding interpolated/extrapolated value \f$y(x)\f$, which can be computed as,
    !>                                      \f{equation*}{
    !>                                           y(x) = \frac {x - x_{i+1}} {x_i - x_{i+1}} y_i + \frac {x - x_{i}} {x_{i+1} - x_{i}} y_{i+1} ~.
    !>                                      \f}
    !>                              <li>    The scalar constant [monopol](@ref pm_polation::monopol) or a scalar object of type [monopol_type](@ref pm_polation::monopol_type)
    !>                                      implying the use of a single **polynomial interpolation/extrapolation** of highest degree `size(crdx) - 1` possible to all pairs of `(crdx, func)` for computing the output `extrap`.<br>
    !>                                      The Neville algorithm is used to approximate the polynomial interpolation/extrapolation.<br>
    !>                          </ol>
    !>  \param[in]  crdx    :   The input `contiguous` vector of<br>
    !>                          <ol>
    !>                              <li>    type `real` of kind \RKALL,
    !>                          </ol>
    !>                          containing the set of abscissa in **strictly ascending order** through which the constructed polynomial must pass.<br>
    !>                          If needed, the input values for `crdx` and `sortedY` can be sorted in ascending order by calling either
    !>                          [setSorted()](@ref pm_arraySort::setSorted) or [setSorted()](@ref pm_arraySort::setSorted).
    !>  \param[in]  func    :   The input `contiguous` vector of the same type, kind, and size as `crdx`,
    !>                          containing the set of function values corresponding to the input abscissa `crdx` through which the constructed polynomial must pass.
    !>  \param[in]  queryx  :   The input scalar or vector of the same type and kind as `crdx`, containing the abscissa (x-value) of the queryx point.
    !>
    !>  \return
    !>
    !>  `extrap`            :   The output object of the same type, kind, rank, and shape as `queryx`,
    !>                          containing the (approximate) interpolated/extrapolated y-value(s) corresponding to the `queryx` point(s).
    !>
    !>  \interface{getExtrap}
    !>  \code{.F90}
    !>
    !>      use pm_polation, only: getExtrap
    !>
    !>      extrap = getExtrap(method, crdx(1:nsam), func(1:nsam), queryx)
    !>      extrap(1:nque) = getExtrap(method, crdx(1:nsam), func(1:nsam), queryx(1:nque))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(crdx) == size(func)` must hold for the corresponding input arguments.<br>
    !>  The condition `isAscendingAll(crdx) .or. same_type_as(method, monopol_type())` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getExtrap](@ref pm_polation::getExtrap)<br>
    !>  [setExtrap](@ref pm_polation::setExtrap)<br>
    !>  [getInterp](@ref pm_polation::getInterp)<br>
    !>  [setInterp](@ref pm_polation::setInterp)<br>
    !>  [pm_sampleQuan](@ref pm_sampleQuan)<br>
    !>  [pm_arraySort](@ref pm_arraySort)<br>
    !>  [pm_quadRomb](@ref pm_quadRomb)<br>
    !>
    !>  \example{getExtrap}
    !>  \include{lineno} example/pm_polation/getExtrap/main.F90
    !>  \compilef{getExtrap}
    !>  \output{getExtrap}
    !>  \include{lineno} example/pm_polation/getExtrap/main.out.F90
    !>  \postproc{getExtrap}
    !>  \include{lineno} example/pm_polation/getExtrap/main.py
    !>  \vis{getExtrap}
    !>  \image html pm_polation/getExtrap/getExtrap.neimean.extrap.png width=700
    !>  \image html pm_polation/getExtrap/getExtrap.neinear.extrap.png width=700
    !>  \image html pm_polation/getExtrap/getExtrap.neinext.extrap.png width=700
    !>  \image html pm_polation/getExtrap/getExtrap.neiprev.extrap.png width=700
    !>  \image html pm_polation/getExtrap/getExtrap.piwilin.extrap.png width=700
    !>  \image html pm_polation/getExtrap/getExtrap.monopol.extrap.png width=700
    !>  \image html pm_polation/getExtrap/getExtrap.rungeEffect.extrap.png width=700
    !>
    !>  \test
    !>  [test_pm_polation](@ref test_pm_polation)
    !>
    !>  \final{getExtrap}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! monopol MNPLD

    interface getExtrap

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getExtrapMNPLD_ND1_QD0_RK5(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapMNPLD_ND1_QD0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(monopol_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getExtrapMNPLD_ND1_QD0_RK4(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapMNPLD_ND1_QD0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(monopol_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getExtrapMNPLD_ND1_QD0_RK3(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapMNPLD_ND1_QD0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(monopol_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getExtrapMNPLD_ND1_QD0_RK2(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapMNPLD_ND1_QD0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(monopol_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getExtrapMNPLD_ND1_QD0_RK1(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapMNPLD_ND1_QD0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(monopol_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getExtrapMNPLD_ND1_QD1_RK5(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapMNPLD_ND1_QD1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(monopol_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getExtrapMNPLD_ND1_QD1_RK4(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapMNPLD_ND1_QD1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(monopol_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getExtrapMNPLD_ND1_QD1_RK3(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapMNPLD_ND1_QD1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(monopol_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getExtrapMNPLD_ND1_QD1_RK2(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapMNPLD_ND1_QD1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(monopol_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getExtrapMNPLD_ND1_QD1_RK1(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapMNPLD_ND1_QD1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(monopol_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! piwilin

    interface getExtrap

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getExtrapPWLN_ND1_QD0_RK5(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapPWLN_ND1_QD0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getExtrapPWLN_ND1_QD0_RK4(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapPWLN_ND1_QD0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getExtrapPWLN_ND1_QD0_RK3(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapPWLN_ND1_QD0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getExtrapPWLN_ND1_QD0_RK2(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapPWLN_ND1_QD0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getExtrapPWLN_ND1_QD0_RK1(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapPWLN_ND1_QD0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getExtrapPWLN_ND1_QD1_RK5(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapPWLN_ND1_QD1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getExtrapPWLN_ND1_QD1_RK4(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapPWLN_ND1_QD1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getExtrapPWLN_ND1_QD1_RK3(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapPWLN_ND1_QD1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getExtrapPWLN_ND1_QD1_RK2(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapPWLN_ND1_QD1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getExtrapPWLN_ND1_QD1_RK1(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapPWLN_ND1_QD1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(piwilin_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! neimean

    interface getExtrap

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getExtrapMEAN_ND1_QD0_RK5(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapMEAN_ND1_QD0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getExtrapMEAN_ND1_QD0_RK4(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapMEAN_ND1_QD0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getExtrapMEAN_ND1_QD0_RK3(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapMEAN_ND1_QD0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getExtrapMEAN_ND1_QD0_RK2(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapMEAN_ND1_QD0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getExtrapMEAN_ND1_QD0_RK1(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapMEAN_ND1_QD0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getExtrapMEAN_ND1_QD1_RK5(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapMEAN_ND1_QD1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getExtrapMEAN_ND1_QD1_RK4(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapMEAN_ND1_QD1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getExtrapMEAN_ND1_QD1_RK3(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapMEAN_ND1_QD1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getExtrapMEAN_ND1_QD1_RK2(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapMEAN_ND1_QD1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getExtrapMEAN_ND1_QD1_RK1(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapMEAN_ND1_QD1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(neimean_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! neinear

    interface getExtrap

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getExtrapNEAR_ND1_QD0_RK5(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapNEAR_ND1_QD0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getExtrapNEAR_ND1_QD0_RK4(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapNEAR_ND1_QD0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getExtrapNEAR_ND1_QD0_RK3(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapNEAR_ND1_QD0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getExtrapNEAR_ND1_QD0_RK2(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapNEAR_ND1_QD0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getExtrapNEAR_ND1_QD0_RK1(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapNEAR_ND1_QD0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getExtrapNEAR_ND1_QD1_RK5(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapNEAR_ND1_QD1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getExtrapNEAR_ND1_QD1_RK4(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapNEAR_ND1_QD1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getExtrapNEAR_ND1_QD1_RK3(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapNEAR_ND1_QD1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getExtrapNEAR_ND1_QD1_RK2(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapNEAR_ND1_QD1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getExtrapNEAR_ND1_QD1_RK1(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapNEAR_ND1_QD1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(neinear_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! neinext

    interface getExtrap

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getExtrapNEXT_ND1_QD0_RK5(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapNEXT_ND1_QD0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getExtrapNEXT_ND1_QD0_RK4(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapNEXT_ND1_QD0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getExtrapNEXT_ND1_QD0_RK3(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapNEXT_ND1_QD0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getExtrapNEXT_ND1_QD0_RK2(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapNEXT_ND1_QD0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getExtrapNEXT_ND1_QD0_RK1(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapNEXT_ND1_QD0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getExtrapNEXT_ND1_QD1_RK5(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapNEXT_ND1_QD1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getExtrapNEXT_ND1_QD1_RK4(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapNEXT_ND1_QD1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getExtrapNEXT_ND1_QD1_RK3(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapNEXT_ND1_QD1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getExtrapNEXT_ND1_QD1_RK2(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapNEXT_ND1_QD1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getExtrapNEXT_ND1_QD1_RK1(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapNEXT_ND1_QD1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(neinext_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! neiprev

    interface getExtrap

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getExtrapPREV_ND1_QD0_RK5(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapPREV_ND1_QD0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getExtrapPREV_ND1_QD0_RK4(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapPREV_ND1_QD0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getExtrapPREV_ND1_QD0_RK3(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapPREV_ND1_QD0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getExtrapPREV_ND1_QD0_RK2(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapPREV_ND1_QD0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getExtrapPREV_ND1_QD0_RK1(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapPREV_ND1_QD0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)                                           :: extrap
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getExtrapPREV_ND1_QD1_RK5(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapPREV_ND1_QD1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    PURE module function getExtrapPREV_ND1_QD1_RK4(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapPREV_ND1_QD1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    PURE module function getExtrapPREV_ND1_QD1_RK3(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapPREV_ND1_QD1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    PURE module function getExtrapPREV_ND1_QD1_RK2(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapPREV_ND1_QD1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    PURE module function getExtrapPREV_ND1_QD1_RK1(method, crdx, func, queryx) result(extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getExtrapPREV_ND1_QD1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)                                           :: extrap(size(queryx, 1, IK))
        type(neiprev_type)  , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the approximate **polynomial interpolation/extrapolation** value of the input specified point `x` for the specified `method`.
    !>
    !>  \details
    !>  For polynomial interpolation/extrapolation, the computation relies on the **Neville algorithm**.<br>
    !>  The extrapolation is done as if the out-of-bound query point is within the boundary (first or last)
    !>  interpolation segment specified by the input `(crdx, func)` pairs of values.<br>
    !>
    !>  \param[in]  method  :   The input scalar constant that can be,
    !>                          <ol>
    !>                              <li>    The scalar constant [neimean](@ref pm_polation::neimean) or a scalar object of type [neimean_type](@ref pm_polation::neimean_type)
    !>                                      implying the use of the average of the `func` values of the two nearest neighbors of the input `queryx` smaller and larger than it as the output `extrap`.<br>
    !>                              <li>    The scalar constant [neinear](@ref pm_polation::neinear) or a scalar object of type [neinear_type](@ref pm_polation::neinear_type)
    !>                                      implying the use of the average of the `func` value of the `neinear` nearest neighbor of the input `queryx` as the output `extrap`.<br>
    !>                                      Note that the nearest neighbor in this case is measured by actual Euclidean distances of neighbors to the input `queryx`.<br>
    !>                              <li>    The scalar constant [neiprev](@ref pm_polation::neiprev) or a scalar object of type [neiprev_type](@ref pm_polation::neiprev_type)
    !>                                      implying the use of the `func` value of the largest abscissa in the input `crdx` smaller than the input `queryx` as the output `extrap`.<br>
    !>                              <li>    The scalar constant [neinext](@ref pm_polation::neinext) or a scalar object of type [neinext_type](@ref pm_polation::neinext_type)
    !>                                      implying the use of the `func` value of the smallest abscissa in the input `crdx` larger than the input `queryx` as the output `extrap`.<br>
    !>                              <li>    The scalar constant [piwilin](@ref pm_polation::piwilin) or a scalar object of type [piwilin_type](@ref pm_polation::piwilin_type)
    !>                                      implying the use of the **linear interpolation/extrapolation** of the `func` values of the two `crdx` points that bracket `queryx` as the output `extrap`.<br>
    !>                                      The linear interpolation/extrapolation implemented in this constructor is based on the Lagrange classical formula for linear interpolation/extrapolation.<br>
    !>                                      Suppose an input query point \f$x\f$ falls between two nodes \f$x_i\f$ and \f$x_{i+1}\f$ with the corresponding function values
    !>                                      \f$y_i\f$ and \f$y_{i+1}\f$ and we wish to estimate the corresponding interpolated/extrapolated value \f$y(x)\f$, which can be computed as,
    !>                                      \f{equation*}{
    !>                                           y(x) = \frac {x - x_{i+1}} {x_i - x_{i+1}} y_i + \frac {x - x_{i}} {x_{i+1} - x_{i}} y_{i+1} ~.
    !>                                      \f}
    !>                              <li>    The scalar constant [monopol](@ref pm_polation::monopol) or a scalar object of type [monopol_type](@ref pm_polation::monopol_type)
    !>                                      implying the use of a single **polynomial interpolation/extrapolation** of highest degree `size(crdx) - 1` possible to all pairs of `(crdx, func)` for computing the output `extrap`.<br>
    !>                                      The Neville algorithm is used to approximate the polynomial interpolation/extrapolation.<br>
    !>                          </ol>
    !>  \param[in]  crdx    :   The input `contiguous` vector of<br>
    !>                          <ol>
    !>                              <li>    type `real` of kind \RKALL,
    !>                          </ol>
    !>                          containing the set of abscissa in **strictly ascending order** through which the constructed polynomial must pass.<br>
    !>                          If needed, the input values for `crdx` and `sortedY` can be sorted in ascending order by calling either
    !>                          [setSorted()](@ref pm_arraySort::setSorted) or [setSorted()](@ref pm_arraySort::setSorted).
    !>  \param[in]  func    :   The input `contiguous` vector of the same type, kind, and size as `crdx`,
    !>                          containing the set of function values corresponding to the input abscissa `crdx` through which the constructed polynomial must pass.
    !>  \param[in]  queryx  :   The input scalar or vector of the same type and kind as `crdx`, containing the abscissa (x-value) of the queryx point.
    !>  \param[out] extrap  :   The output object of the same type, kind, rank, and shape as `queryx`,
    !>                          containing the (approximate) interpolated/extrapolated y-value(s) corresponding to the `queryx` point(s).
    !>  \param[out] relerr  :   The output scalar of the same type and kind as `crdx`, containing the estimated error in the output `extrap`.
    !>                          (**optional**. It can be present **if and only** the input argument `method` is set to
    !>                          [monopol](@ref pm_polation::monopol) or a scalar object of type [monopol_type](@ref pm_polation::monopol_type) .)
    !>
    !>  \interface{setExtrap}
    !>  \code{.F90}
    !>
    !>      use pm_polation, only: setExtrap
    !>
    !>      call setExtrap(method, crdx(1:nsam), func(1:nsam), queryx, extrap)
    !>      call setExtrap(method, crdx(1:nsam), func(1:nsam), queryx(1:nque), extrap(1:nque))
    !>
    !>      call setExtrap(method, crdx(1:nsam), func(1:nsam), queryx, extrap, relerr) ! method = monopol
    !>      call setExtrap(method, crdx(1:nsam), func(1:nsam), queryx(1:nque), extrap(1:nque), relerr(1:nque)) ! method = monopol
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(crdx) == size(func)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(shape(queryx) == shape(extrap))` must hold for the corresponding input arguments.<br>
    !>  The condition `all(shape(queryx) == shape(relerr))` must hold for the corresponding input arguments.<br>
    !>  The condition `isAscendingAll(crdx) .or. same_type_as(method, monopol_type())` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getExtrap](@ref pm_polation::getExtrap)<br>
    !>  [setExtrap](@ref pm_polation::setExtrap)<br>
    !>  [getInterp](@ref pm_polation::getInterp)<br>
    !>  [setInterp](@ref pm_polation::setInterp)<br>
    !>  [pm_sampleQuan](@ref pm_sampleQuan)<br>
    !>  [pm_arraySort](@ref pm_arraySort)<br>
    !>  [pm_quadRomb](@ref pm_quadRomb)<br>
    !>
    !>  \example{setExtrap}
    !>  \include{lineno} example/pm_polation/setExtrap/main.F90
    !>  \compilef{setExtrap}
    !>  \output{setExtrap}
    !>  \include{lineno} example/pm_polation/setExtrap/main.out.F90
    !>  \postproc{setExtrap}
    !>  \include{lineno} example/pm_polation/setExtrap/main.py
    !>  \vis{setExtrap}
    !>  \image html pm_polation/setExtrap/setExtrap.neimean.extrap.png width=700
    !>  \image html pm_polation/setExtrap/setExtrap.neinear.extrap.png width=700
    !>  \image html pm_polation/setExtrap/setExtrap.neinext.extrap.png width=700
    !>  \image html pm_polation/setExtrap/setExtrap.neiprev.extrap.png width=700
    !>  \image html pm_polation/setExtrap/setExtrap.piwilin.extrap.png width=700
    !>  \image html pm_polation/setExtrap/setExtrap.monopol.extrap.png width=700
    !>  \image html pm_polation/setExtrap/setExtrap.rungeEffect.extrap.png width=700
    !>
    !>  \test
    !>  [test_pm_polation](@ref test_pm_polation)
    !>
    !>  \final{setExtrap}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! monopol MNPLD

    interface setExtrap

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setExtrapMNPLD_ND1_QD0_RK5(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMNPLD_ND1_QD0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setExtrapMNPLD_ND1_QD0_RK4(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMNPLD_ND1_QD0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setExtrapMNPLD_ND1_QD0_RK3(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMNPLD_ND1_QD0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setExtrapMNPLD_ND1_QD0_RK2(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMNPLD_ND1_QD0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setExtrapMNPLD_ND1_QD0_RK1(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMNPLD_ND1_QD0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setExtrapMNPLD_ND1_QD1_RK5(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMNPLD_ND1_QD1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setExtrapMNPLD_ND1_QD1_RK4(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMNPLD_ND1_QD1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setExtrapMNPLD_ND1_QD1_RK3(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMNPLD_ND1_QD1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setExtrapMNPLD_ND1_QD1_RK2(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMNPLD_ND1_QD1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setExtrapMNPLD_ND1_QD1_RK1(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMNPLD_ND1_QD1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! monopol MNPLE

    interface setExtrap

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setExtrapMNPLE_ND1_QD0_RK5(method, crdx, func, queryx, extrap, relerr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMNPLE_ND1_QD0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap, relerr
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setExtrapMNPLE_ND1_QD0_RK4(method, crdx, func, queryx, extrap, relerr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMNPLE_ND1_QD0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap, relerr
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setExtrapMNPLE_ND1_QD0_RK3(method, crdx, func, queryx, extrap, relerr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMNPLE_ND1_QD0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap, relerr
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setExtrapMNPLE_ND1_QD0_RK2(method, crdx, func, queryx, extrap, relerr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMNPLE_ND1_QD0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap, relerr
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setExtrapMNPLE_ND1_QD0_RK1(method, crdx, func, queryx, extrap, relerr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMNPLE_ND1_QD0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap, relerr
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setExtrapMNPLE_ND1_QD1_RK5(method, crdx, func, queryx, extrap, relerr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMNPLE_ND1_QD1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:), relerr(:)
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setExtrapMNPLE_ND1_QD1_RK4(method, crdx, func, queryx, extrap, relerr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMNPLE_ND1_QD1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:), relerr(:)
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setExtrapMNPLE_ND1_QD1_RK3(method, crdx, func, queryx, extrap, relerr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMNPLE_ND1_QD1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:), relerr(:)
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setExtrapMNPLE_ND1_QD1_RK2(method, crdx, func, queryx, extrap, relerr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMNPLE_ND1_QD1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:), relerr(:)
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setExtrapMNPLE_ND1_QD1_RK1(method, crdx, func, queryx, extrap, relerr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMNPLE_ND1_QD1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:), relerr(:)
        type(monopol_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! piwilin

    interface setExtrap

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setExtrapPWLN_ND1_QD0_RK5(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapPWLN_ND1_QD0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(piwilin_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setExtrapPWLN_ND1_QD0_RK4(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapPWLN_ND1_QD0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(piwilin_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setExtrapPWLN_ND1_QD0_RK3(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapPWLN_ND1_QD0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(piwilin_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setExtrapPWLN_ND1_QD0_RK2(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapPWLN_ND1_QD0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(piwilin_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setExtrapPWLN_ND1_QD0_RK1(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapPWLN_ND1_QD0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(piwilin_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setExtrapPWLN_ND1_QD1_RK5(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapPWLN_ND1_QD1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(piwilin_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setExtrapPWLN_ND1_QD1_RK4(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapPWLN_ND1_QD1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(piwilin_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setExtrapPWLN_ND1_QD1_RK3(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapPWLN_ND1_QD1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(piwilin_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setExtrapPWLN_ND1_QD1_RK2(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapPWLN_ND1_QD1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(piwilin_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setExtrapPWLN_ND1_QD1_RK1(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapPWLN_ND1_QD1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(piwilin_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! neimean

    interface setExtrap

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setExtrapMEAN_ND1_QD0_RK5(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMEAN_ND1_QD0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(neimean_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setExtrapMEAN_ND1_QD0_RK4(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMEAN_ND1_QD0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(neimean_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setExtrapMEAN_ND1_QD0_RK3(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMEAN_ND1_QD0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(neimean_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setExtrapMEAN_ND1_QD0_RK2(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMEAN_ND1_QD0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(neimean_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setExtrapMEAN_ND1_QD0_RK1(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMEAN_ND1_QD0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(neimean_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setExtrapMEAN_ND1_QD1_RK5(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMEAN_ND1_QD1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(neimean_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setExtrapMEAN_ND1_QD1_RK4(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMEAN_ND1_QD1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(neimean_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setExtrapMEAN_ND1_QD1_RK3(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMEAN_ND1_QD1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(neimean_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setExtrapMEAN_ND1_QD1_RK2(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMEAN_ND1_QD1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(neimean_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setExtrapMEAN_ND1_QD1_RK1(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapMEAN_ND1_QD1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(neimean_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! neinear

    interface setExtrap

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setExtrapNEAR_ND1_QD0_RK5(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapNEAR_ND1_QD0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(neinear_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setExtrapNEAR_ND1_QD0_RK4(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapNEAR_ND1_QD0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(neinear_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setExtrapNEAR_ND1_QD0_RK3(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapNEAR_ND1_QD0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(neinear_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setExtrapNEAR_ND1_QD0_RK2(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapNEAR_ND1_QD0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(neinear_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setExtrapNEAR_ND1_QD0_RK1(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapNEAR_ND1_QD0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(neinear_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setExtrapNEAR_ND1_QD1_RK5(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapNEAR_ND1_QD1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(neinear_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setExtrapNEAR_ND1_QD1_RK4(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapNEAR_ND1_QD1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(neinear_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setExtrapNEAR_ND1_QD1_RK3(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapNEAR_ND1_QD1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(neinear_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setExtrapNEAR_ND1_QD1_RK2(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapNEAR_ND1_QD1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(neinear_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setExtrapNEAR_ND1_QD1_RK1(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapNEAR_ND1_QD1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(neinear_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! neinext

    interface setExtrap

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setExtrapNEXT_ND1_QD0_RK5(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapNEXT_ND1_QD0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(neinext_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setExtrapNEXT_ND1_QD0_RK4(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapNEXT_ND1_QD0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(neinext_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setExtrapNEXT_ND1_QD0_RK3(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapNEXT_ND1_QD0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(neinext_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setExtrapNEXT_ND1_QD0_RK2(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapNEXT_ND1_QD0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(neinext_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setExtrapNEXT_ND1_QD0_RK1(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapNEXT_ND1_QD0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(neinext_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setExtrapNEXT_ND1_QD1_RK5(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapNEXT_ND1_QD1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(neinext_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setExtrapNEXT_ND1_QD1_RK4(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapNEXT_ND1_QD1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(neinext_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setExtrapNEXT_ND1_QD1_RK3(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapNEXT_ND1_QD1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(neinext_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setExtrapNEXT_ND1_QD1_RK2(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapNEXT_ND1_QD1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(neinext_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setExtrapNEXT_ND1_QD1_RK1(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapNEXT_ND1_QD1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(neinext_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! neiprev

    interface setExtrap

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setExtrapPREV_ND1_QD0_RK5(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapPREV_ND1_QD0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(neiprev_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setExtrapPREV_ND1_QD0_RK4(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapPREV_ND1_QD0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(neiprev_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setExtrapPREV_ND1_QD0_RK3(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapPREV_ND1_QD0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(neiprev_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setExtrapPREV_ND1_QD0_RK2(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapPREV_ND1_QD0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(neiprev_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setExtrapPREV_ND1_QD0_RK1(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapPREV_ND1_QD0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)                    :: queryx
        real(RKG)           , intent(out)                   :: extrap
        type(neiprev_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setExtrapPREV_ND1_QD1_RK5(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapPREV_ND1_QD1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(neiprev_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setExtrapPREV_ND1_QD1_RK4(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapPREV_ND1_QD1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(neiprev_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setExtrapPREV_ND1_QD1_RK3(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapPREV_ND1_QD1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(neiprev_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setExtrapPREV_ND1_QD1_RK2(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapPREV_ND1_QD1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(neiprev_type)  , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setExtrapPREV_ND1_QD1_RK1(method, crdx, func, queryx, extrap)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setExtrapPREV_ND1_QD1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: crdx(:), func(:)
        real(RKG)           , intent(in)    , contiguous    :: queryx(:)
        real(RKG)           , intent(out)   , contiguous    :: extrap(:)
        type(neiprev_type)  , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_polation