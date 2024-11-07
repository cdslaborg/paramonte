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
!>  This module contains procedures, generic interfaces, and types for numerical optimizations of mathematical functions.
!>
!>  \details
!>  The following methods are currently included in this module:
!>  <ol>
!>      <li>    [setBracketMin](@ref pm_optimization::setBracketMin), creates an interval triplet guaranteed to contain the minimum of convex function.
!>      <li>    [setBracketMax](@ref pm_optimization::setBracketMax), creates an interval triplet guaranteed to contain the maximum of convex function.
!>      <li>    [isBracketMin](@ref pm_optimization::isBracketMin), tests if an interval triplet contains the minimum of convex function.
!>      <li>    [isBracketMax](@ref pm_optimization::isBracketMax), tests if an interval triplet contains the maximum of convex function.
!>      <li>    [getMinBrent](@ref pm_optimization::getMinBrent) and [setMinBrent](@ref pm_optimization::setMinBrent) perform unconstrained derivative-free minimization of scalar **univariate** functions.<br>
!>              **The Brent method** is a hybrid iterative minimization algorithm combining the bisection method, the secant method and inverse quadratic interpolation.<br>
!>              It has the reliability of bisection or Golden Section methods but it can be as quick as some of the less-reliable methods.<br>
!>              The algorithm tries to use the potentially fast-converging secant method or inverse quadratic interpolation if possible,
!>              but it falls back to the more robust bisection method if necessary.<br>
!>              The Brent method is due to [Richard Brent](https://en.wikipedia.org/wiki/Richard_P._Brent) and builds on an earlier algorithm by [Theodorus Dekker](https://en.wikipedia.org/wiki/Theodorus_Dekker).<br>
!>              Consequently, the method is also known as the Brent–Dekker method.<br>
!>      <li>    [isFailedMinPowell](@ref pm_optimization::isFailedMinPowell) and [setMinPowell](@ref pm_optimization::setMinPowell) perform unconstrained derivative-free minimization of **multivariate** functions.<br>
!>              **The Powell method**, strictly the **Powell conjugate direction method**, is an algorithm proposed by
!>              [Michael J. D. Powell](https://en.wikipedia.org/wiki/Michael_J._D._Powell) for finding a local minimum of a function.<br>
!>              The function need not be differentiable, and no derivatives are taken.<br>
!>              The function must be a real-valued function of a fixed number of real-valued inputs.<br>
!>              The caller passes in the initial point.<br>
!>              The caller also passes in a set of initial search vectors.<br>
!>              Typically \f$N\f$ search vectors \f$\{s_{1}, \dots, s_{N}\}\f$ are passed in, which are simply the normal vectors aligned to each axis.<br>
!>              The method minimizes the function by a bi-directional search along each search vector, in turn.<br>
!>              The bi-directional line search along each search vector can be done by the Golden-section search or the Brent method.<br>
!>              Let the minima found during each bi-directional line search be \f$\{x_{0}+\alpha_{1}s_{1},{x}_{0} + \sum_{i=1}^{2}\alpha_{i}{s}_{i},\dots ,{x}_{0}+\sum_{i=1}^{N}\alpha_{i}{s}_{i}\}\f$,
!>              where \f$x_0\f$ is the initial starting point and \f$\alpha_i\f$ is the scalar determined during bi-directional search along \f$s_i\f$.<br>
!>              The new position \f$x_{1}\f$ can then be expressed as a linear combination of the search vectors, i.e., \f$x_{1} = x_{0} + \sum_{i=1}^{N} \alpha_{i}s_{i}\f$.<br>
!>              The new displacement vector \f$\sum_{i = 1}^{N} \alpha_{i}s_{i}\f$ becomes a new search vector, and is added to the end of the search vector list.<br>
!>              Meanwhile, the search vector which contributed most to the new direction, i.e. the one which was most successful \f$i_{d} = \arg \max_{i = 1}^{N} | \alpha_{i} | \| s_{i} \|\f$, is deleted from the search vector list.<br>
!>              The new set of \f$N\f$ search vectors is \f$\{s_{1}, \dots , s_{i_{d} - 1}, s_{i_{d} + 1}, \dots , s_{N}, \sum_{i=1}^{N} \alpha_{i}s_{i}\}\f$.<br>
!>              The algorithm iterates an arbitrary number of times until no significant improvement is made.<br>
!>              The method is useful for calculating the local minimum of a continuous but complex function,
!>              especially one without an underlying mathematical definition, because it is not necessary to take derivatives.<br>
!>              The basic algorithm is simple and the complexity is in the linear searches along the search vectors, which can be achieved via the Brent method.<br>
!       <li>    [getMinSimplex](@ref pm_optimization::getMinSimplex) and [setMinSimplex](@ref pm_optimization::setMinSimplex) perform unconstrained derivative-free minimization of **multivariate** functions.<br>
!               The **Nelder–Mead method** (also **downhill simplex method**, **amoeba method**, or **polytope method**)
!               is a numerical method used to find the minimum or maximum of an objective function in a **multidimensional** space.<br>
!               It is a direct search method (based on function comparison) and is often applied to **nonlinear optimization problems** for which **derivatives are unknown**.<br>
!               The Nelder–Mead technique is a **heuristic search method** that can converge to **non-stationary points**.<br>
!               The Nelder–Mead technique was proposed by [John Nelder](https://en.wikipedia.org/wiki/John_Nelder) and [Roger Mead](https://en.wikipedia.org/wiki/Roger_Mead) in 1965.<br>
!               According to Numerical Recipes by Press et al. (1992):<br>
!               >  The Nelder–Mead method crawls downhill in a straightforward fashion that makes almost no special assumptions about the target function.<br>
!               >  This can be extremely slow, but it can also, in some cases, be extremely robust.<br>
!               >  The storage requirement for this method is of order \f$\ms{ndim}^2\f$, where \f$\ms{ndim}\f$ is the number of dimensions of the function support.<br>
!               >  As such, the Nedler-Mead Simplex method is most useful when the minimization tasks are minimal and isolated.<br>
!>  </ol>
!>
!>  \test
!>  [test_pm_optimization](@ref test_pm_optimization)
!>  
!>  \todo
!>  \phigh
!>  The Nedler-Mead algorithm must be implemented.<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_optimization

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_optimization"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if and only if a concave quadratic curve can fit the specified input triple [xmin, xlow, xupp]`
    !>  and the function value at the middle point `xmin` is larger than both boundary point function values.<br>
    !>
    !>  \details
    !>  This function is useful for testing whether the initial bracket intervals
    !>  to various optimization routines enclose the desired function extremum.<br>
    !>
    !>  Specifically, this generic interface computed and returns the following criterion:
    !>  `bracketed = xlow <= xmax .and. xmax <= xupp .and. flow <= fmax .and. fmax >= fupp`
    !>
    !>  \param[in]      xmax    :   The input scalar of type `real` of kind \RKALL,
    !>                              containing the abscissa corresponding to the maximum function value among the triplet `[xmax, xlow, xupp]`.<br>
    !>  \param[in]      xlow    :   The input scalar of the same type and kind as the input argument `xmax`,
    !>                              containing the lower bound of the interval specified by the triplet.<br>
    !>  \param[in]      xupp    :   The input scalar of the same type and kind as the input argument `xmax`,
    !>                              containing the upper bound of the interval specified by the triplet.<br>
    !>  \param[in]      fmax    :   The input scalar of the same type and kind as the input argument `xmax`,
    !>                              containing the target function value at `xmax`.<br>
    !>  \param[in]      flow    :   The input scalar of the same type and kind as the input argument `xmax`,
    !>                              containing the target function value at `xlow`.<br>
    !>  \param[in]      fupp    :   The input scalar of the same type and kind as the input argument `xmax`,
    !>                              containing the target function value at `xupp`.<br>
    !>
    !>  \interface{isBracketMax}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_optimization, only: isBracketMax
    !>      logical(LK) :: bracketed
    !>
    !>      bracketed = isBracketMax(xmax, xlow, xupp, fmax, flow, fupp)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [isBracketMax](@ref pm_optimization::isBracketMax)<br>
    !>  [isBracketMin](@ref pm_optimization::isBracketMin)<br>
    !>  [setBracketMax](@ref pm_optimization::setBracketMax)<br>
    !>  [setBracketMin](@ref pm_optimization::setBracketMin)<br>
    !>  [getMinBrent](@ref pm_optimization::getMinBrent)<br>
    !>  [setMinBrent](@ref pm_optimization::setMinBrent)<br>
    !>
    !>  \example{isBracketMax}
    !>  \include{lineno} example/pm_optimization/isBracketMax/main.F90
    !>  \compilef{isBracketMax}
    !>  \output{isBracketMax}
    !>  \include{lineno} example/pm_optimization/isBracketMax/main.out.F90
    !>
    !>  \test
    !>  [test_pm_optimization](@ref test_pm_optimization)
    !>
    !>  \final{isBracketMax}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isBracketMax

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function isBracketMax_RK5(xmax, xlow, xupp, fmax, flow, fupp) result(bracketed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isBracketMax_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)                :: xmax, xlow, xupp, fmax, flow, fupp
        logical(LK)                                     :: bracketed
    end function
#endif

#if RK4_ENABLED
    pure elemental module function isBracketMax_RK4(xmax, xlow, xupp, fmax, flow, fupp) result(bracketed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isBracketMax_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)                :: xmax, xlow, xupp, fmax, flow, fupp
        logical(LK)                                     :: bracketed
    end function
#endif

#if RK3_ENABLED
    pure elemental module function isBracketMax_RK3(xmax, xlow, xupp, fmax, flow, fupp) result(bracketed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isBracketMax_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)                :: xmax, xlow, xupp, fmax, flow, fupp
        logical(LK)                                     :: bracketed
    end function
#endif

#if RK2_ENABLED
    pure elemental module function isBracketMax_RK2(xmax, xlow, xupp, fmax, flow, fupp) result(bracketed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isBracketMax_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)                :: xmax, xlow, xupp, fmax, flow, fupp
        logical(LK)                                     :: bracketed
    end function
#endif

#if RK1_ENABLED
    pure elemental module function isBracketMax_RK1(xmax, xlow, xupp, fmax, flow, fupp) result(bracketed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isBracketMax_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)                :: xmax, xlow, xupp, fmax, flow, fupp
        logical(LK)                                     :: bracketed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if and only if a convex quadratic curve can fit the specified input triple `[xmin, xlow, xupp]`
    !>  and the function value at the middle point `xmin` is smaller than both boundary point function values.<br>
    !>
    !>  \details
    !>  This function is useful for testing whether the initial bracket intervals
    !>  to various optimization routines enclose the desired function extremum.<br>
    !>
    !>  Specifically, this generic interface computed and returns the following criterion:
    !>  `bracketed = xlow <= xmin .and. xmin <= xupp .and. flow >= fmin .or. fmin <= fupp`
    !>
    !>  \param[in]      xmin    :   The input scalar of type `real` of kind \RKALL,
    !>                              containing the abscissa corresponding to the minimum function value among the triplet `[xmin, xlow, xupp]`.<br>
    !>  \param[in]      xlow    :   The input scalar of the same type and kind as the input argument `xmin`,
    !>                              containing the lower bound of the interval specified by the triplet.<br>
    !>  \param[in]      xupp    :   The input scalar of the same type and kind as the input argument `xmin`,
    !>                              containing the upper bound of the interval specified by the triplet.<br>
    !>  \param[in]      fmin    :   The input scalar of the same type and kind as the input argument `xmin`,
    !>                              containing the target function value at `xmin`.<br>
    !>  \param[in]      flow    :   The input scalar of the same type and kind as the input argument `xmin`,
    !>                              containing the target function value at `xlow`.<br>
    !>  \param[in]      fupp    :   The input scalar of the same type and kind as the input argument `xmin`,
    !>                              containing the target function value at `xupp`.<br>
    !>
    !>  \interface{isBracketMin}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_optimization, only: isBracketMin
    !>      logical(LK) :: bracketed
    !>
    !>      bracketed = isBracketMin(xmin, xlow, xupp, fmin, flow, fupp)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [isBracketMax](@ref pm_optimization::isBracketMax)<br>
    !>  [isBracketMin](@ref pm_optimization::isBracketMin)<br>
    !>  [setBracketMax](@ref pm_optimization::setBracketMax)<br>
    !>  [setBracketMin](@ref pm_optimization::setBracketMin)<br>
    !>  [getMinBrent](@ref pm_optimization::getMinBrent)<br>
    !>  [setMinBrent](@ref pm_optimization::setMinBrent)<br>
    !>
    !>  \example{isBracketMin}
    !>  \include{lineno} example/pm_optimization/isBracketMin/main.F90
    !>  \compilef{isBracketMin}
    !>  \output{isBracketMin}
    !>  \include{lineno} example/pm_optimization/isBracketMin/main.out.F90
    !>
    !>  \test
    !>  [test_pm_optimization](@ref test_pm_optimization)
    !>
    !>  \final{isBracketMin}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isBracketMin

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function isBracketMin_RK5(xmin, xlow, xupp, fmin, flow, fupp) result(bracketed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isBracketMin_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)                :: xmin, xlow, xupp, fmin, flow, fupp
        logical(LK)                                     :: bracketed
    end function
#endif

#if RK4_ENABLED
    pure elemental module function isBracketMin_RK4(xmin, xlow, xupp, fmin, flow, fupp) result(bracketed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isBracketMin_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)                :: xmin, xlow, xupp, fmin, flow, fupp
        logical(LK)                                     :: bracketed
    end function
#endif

#if RK3_ENABLED
    pure elemental module function isBracketMin_RK3(xmin, xlow, xupp, fmin, flow, fupp) result(bracketed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isBracketMin_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)                :: xmin, xlow, xupp, fmin, flow, fupp
        logical(LK)                                     :: bracketed
    end function
#endif

#if RK2_ENABLED
    pure elemental module function isBracketMin_RK2(xmin, xlow, xupp, fmin, flow, fupp) result(bracketed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isBracketMin_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)                :: xmin, xlow, xupp, fmin, flow, fupp
        logical(LK)                                     :: bracketed
    end function
#endif

#if RK1_ENABLED
    pure elemental module function isBracketMin_RK1(xmin, xlow, xupp, fmin, flow, fupp) result(bracketed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isBracketMin_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)                :: xmin, xlow, xupp, fmin, flow, fupp
        logical(LK)                                     :: bracketed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Refine an initial input interval such that the final returned interval is
    !>  guaranteed to contain the maximum of the user-specified input function.<br>
    !>
    !>  \details
    !>  Given two initial input values `[xlow, xupp]` representing the best-guess interval to hopefully contain the maximum of the input function,
    !>  compute and return the refined interval and an interior point `xmax` within the interval such that the triple `[xlow, xmax, xupp]` are
    !>  guaranteed to pass through a convex parabola.<br>
    !>  In addition, return the function values at the final computed triplet.<br>
    !>  The output refined interval `[xlow, xupp]` is guaranteed to contain a maximum of the user-specified input function.<br>
    !>
    !>  \warning
    !>  The final output `xmax` is **not** the location of a function maximum.<br>
    !>  It merely satisfies the condition `xlow < xmax .and. xmax < xupp`.<br>
    !>
    !>  \param[in]      getFunc :   The scalar function be maximized.
    !>                              <ol>
    !>                                  <li>    On input, it must take a scalar of the same type and kind as the output argument `xmax`.<br>
    !>                                  <li>    On output, it must return a scalar of the same type and kind as the output argument `xmax`,
    !>                                          containing the function value at the specified input scalar point.<br>
    !>                              </ol>
    !>                              The following demonstrates the interface of `getFunc`,
    !>                              \code{.F90}
    !>                                  function getFunc(x) result(func)
    !>                                      real(RKG), intent(in) :: x
    !>                                      real(RKG) :: func
    !>                                  end function
    !>                              \endcode
    !>                              where `RKG` can refer to any `real` type kind parameter \RKALL supported by the library.<br>
    !>  \param[inout]   niter   :   The input/output positive scalar of type `integer` of default kind \IK containing
    !>                              the **maximum** number of allowed iterations in the algorithm in search of the bracket.<br>
    !>                              On output,
    !>                              <ol>
    !>                                  <li>    **If the algorithm succeeds**, `niter` will be set to the actual number of iterations
    !>                                          taken to find the bracket which is by definition, less than or equal to the input value.<br>
    !>                                  <li>    **If the algorithm fails to converge**, it will be a number larger than the input value (by only `1` unit).<br>
    !>                              </ol>
    !>                              The value of `niter` is almost twice the maximum allowed number of calls to the user-specified input function.<br>
    !>                              The reasonable choice can be `1000`.<br>
    !>  \param[out]     xmax    :   The output scalar of type `real` of kind \RKALL,
    !>                              containing the abscissa corresponding to the maximum function value among the triplet `[xlow, xmax, xupp]`.<br>
    !>  \param[inout]   xlow    :   The input/output scalar of the same type and kind as the output argument `xmax`,
    !>                              containing the lower bound of the search interval for the function maximum abscissa.<br>
    !>                              The condition `xlow < xmax .and. getFunc(xmax) < getFunc(xlow)` must hold for the corresponding input arguments.<br>
    !>  \param[inout]   xupp    :   The input/output scalar of the same type and kind as the output argument `xmax`,
    !>                              containing the upper bound of the search interval for the function maximum abscissa.<br>
    !>                              The condition `xmax < xupp .and. getFunc(xmax) < getFunc(xupp)` must hold for the corresponding input arguments.<br>
    !>  \param[out]     fmax    :   The output scalar of the same type and kind as the output argument `xmax`,
    !>                              containing the user-specified function value at `xmax`.<br>
    !>  \param[out]     flow    :   The output scalar of the same type and kind as the output argument `xmax`,
    !>                              containing the user-specified function value at the final refined `xlow`.<br>
    !>                              (**optional**.)
    !>  \param[out]     fupp    :   The output scalar of the same type and kind as the output argument `xmax`,
    !>                              containing the user-specified function value at the final refined `xupp`.<br>
    !>                              (**optional**.)
    !>
    !>  \interface{setBracketMax}
    !>  \code{.F90}
    !>
    !>      use pm_optimization, only: setBracketMax
    !>
    !>      call setBracketMax(getFunc, niter, xmax, xlow, xupp, fmax, flow, fupp)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < niter` must hold for the corresponding input arguments.<br>
    !>  The condition `xlow < xupp` must hold for the corresponding input arguments.<br>
    !>  The input function `getFunc()` must be concave (with negative second derivative) at least in one interval within its domain,
    !>  otherwise the procedure will enter a semi-infinite loop in search of the maximum.<br>
    !>
    !>  \impure
    !>
    !>  \recursive
    !>
    !>  \see
    !>  [getMinBrent](@ref pm_optimization::getMinBrent)<br>
    !>  [setMinBrent](@ref pm_optimization::setMinBrent)<br>
    !>
    !>  \example{setBracketMax}
    !>  \include{lineno} example/pm_optimization/setBracketMax/main.F90
    !>  \compilef{setBracketMax}
    !>  \output{setBracketMax}
    !>  \include{lineno} example/pm_optimization/setBracketMax/main.out.F90
    !>
    !>  \test
    !>  [test_pm_optimization](@ref test_pm_optimization)
    !>
    !>  \final{setBracketMax}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setBracketMax

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module subroutine setBracketMax_RK5(getFunc, niter, xmax, xlow, xupp, fmax, flow, fupp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBracketMax_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                            :: getFunc
        integer(IK)         , intent(inout)             :: niter
        real(RKG)           , intent(inout)             :: xlow, xupp
        real(RKG)           , intent(out)               :: xmax, fmax
        real(RKG)           , intent(out)   , optional  :: flow, fupp
    end subroutine
#endif

#if RK4_ENABLED
    recursive module subroutine setBracketMax_RK4(getFunc, niter, xmax, xlow, xupp, fmax, flow, fupp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBracketMax_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                            :: getFunc
        integer(IK)         , intent(inout)             :: niter
        real(RKG)           , intent(inout)             :: xlow, xupp
        real(RKG)           , intent(out)               :: xmax, fmax
        real(RKG)           , intent(out)   , optional  :: flow, fupp
    end subroutine
#endif

#if RK3_ENABLED
    recursive module subroutine setBracketMax_RK3(getFunc, niter, xmax, xlow, xupp, fmax, flow, fupp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBracketMax_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                            :: getFunc
        integer(IK)         , intent(inout)             :: niter
        real(RKG)           , intent(inout)             :: xlow, xupp
        real(RKG)           , intent(out)               :: xmax, fmax
        real(RKG)           , intent(out)   , optional  :: flow, fupp
    end subroutine
#endif

#if RK2_ENABLED
    recursive module subroutine setBracketMax_RK2(getFunc, niter, xmax, xlow, xupp, fmax, flow, fupp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBracketMax_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                            :: getFunc
        integer(IK)         , intent(inout)             :: niter
        real(RKG)           , intent(inout)             :: xlow, xupp
        real(RKG)           , intent(out)               :: xmax, fmax
        real(RKG)           , intent(out)   , optional  :: flow, fupp
    end subroutine
#endif

#if RK1_ENABLED
    recursive module subroutine setBracketMax_RK1(getFunc, niter, xmax, xlow, xupp, fmax, flow, fupp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBracketMax_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                            :: getFunc
        integer(IK)         , intent(inout)             :: niter
        real(RKG)           , intent(inout)             :: xlow, xupp
        real(RKG)           , intent(out)               :: xmax, fmax
        real(RKG)           , intent(out)   , optional  :: flow, fupp
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Refine an initial input interval such that the final returned interval is
    !>  guaranteed to contain the minimum of the user-specified input function.<br>
    !>
    !>  \details
    !>  Given two initial input values `[xlow, xupp]` representing the best-guess interval to hopefully contain the minimum of the input function,
    !>  compute and return the refined interval and an interior point `xmin` within the interval such that the triple `[xlow, xmin, xupp]` are
    !>  guaranteed to pass through a convex parabola.<br>
    !>  In addition, return the function values at the final computed triplet.<br>
    !>  The output refined interval `[xlow, xupp]` is guaranteed to contain a minimum of the user-specified input function.<br>
    !>
    !>  \warning
    !>  The final output `xmin` is **not** the location of a function minimum.<br>
    !>  It merely satisfies the condition `xlow < xmin .and. xmin < xupp`.<br>
    !>
    !>  \param[in]      getFunc :   The scalar function be minimized.
    !>                              <ol>
    !>                                  <li>    On input, it must take a scalar of the same type and kind as the output argument `xmin`.<br>
    !>                                  <li>    On output, it must return a scalar of the same type and kind as the output argument `xmin`,
    !>                                          containing the function value at the specified input scalar point.<br>
    !>                              </ol>
    !>                              The following demonstrates the interface of `getFunc`,
    !>                              \code{.F90}
    !>                                  function getFunc(x) result(func)
    !>                                      real(RKG), intent(in) :: x
    !>                                      real(RKG) :: func
    !>                                  end function
    !>                              \endcode
    !>                              where `RKG` can refer to any `real` type kind parameter \RKALL supported by the library.<br>
    !>  \param[inout]   niter   :   The input/output positive scalar of type `integer` of default kind \IK containing
    !>                              the **maximum** number of allowed iterations in the algorithm in search of the bracket.<br>
    !>                              On output,
    !>                              <ol>
    !>                                  <li>    **If the algorithm succeeds**, `niter` will be set to the actual number of iterations
    !>                                          taken to find the bracket which is by definition, less than or equal to the input value.<br>
    !>                                  <li>    **If the algorithm fails to converge**, it will be a number larger than the input value (by only `1` unit).<br>
    !>                              </ol>
    !>                              The value of `niter` is almost twice the maximum allowed number of calls to the user-specified input function.<br>
    !>                              The reasonable choice can be `1000`.<br>
    !>  \param[out]     xmin    :   The output scalar of type `real` of kind \RKALL,
    !>                              containing the abscissa corresponding to the minimum function value among the triplet `[xlow, xmin, xupp]`.<br>
    !>  \param[inout]   xlow    :   The input/output scalar of the same type and kind as the output argument `xmin`,
    !>                              containing the lower bound of the search interval for the function minimum abscissa.<br>
    !>                              The condition `xlow < xmin .and. getFunc(xmin) < getFunc(xlow)` must hold for the corresponding input arguments.<br>
    !>  \param[inout]   xupp    :   The input/output scalar of the same type and kind as the output argument `xmin`,
    !>                              containing the upper bound of the search interval for the function minimum abscissa.<br>
    !>                              The condition `xmin < xupp .and. getFunc(xmin) < getFunc(xupp)` must hold for the corresponding input arguments.<br>
    !>  \param[out]     fmin    :   The output scalar of the same type and kind as the output argument `xmin`,
    !>                              containing the user-specified function value at `xmin`.<br>
    !>  \param[out]     flow    :   The output scalar of the same type and kind as the output argument `xmin`,
    !>                              containing the user-specified function value at the final refined `xlow`.<br>
    !>                              (**optional**.)
    !>  \param[out]     fupp    :   The output scalar of the same type and kind as the output argument `xmin`,
    !>                              containing the user-specified function value at the final refined `xupp`.<br>
    !>                              (**optional**.)
    !>
    !>  \interface{setBracketMin}
    !>  \code{.F90}
    !>
    !>      use pm_optimization, only: setBracketMin
    !>
    !>      call setBracketMin(getFunc, niter, xmin, xlow, xupp, fmin, flow, fupp)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < niter` must hold for the corresponding input arguments.<br>
    !>  The condition `xlow < xupp` must hold for the corresponding input arguments.<br>
    !>  The input function `getFunc()` must be convex (with positive second derivative) at least in one interval within its domain,
    !>  otherwise the procedure will enter a semi-infinite loop in search of the minimum.<br>
    !>
    !>  \impure
    !>
    !>  \recursive
    !>
    !>  \see
    !>  [getMinBrent](@ref pm_optimization::getMinBrent)<br>
    !>  [setMinBrent](@ref pm_optimization::setMinBrent)<br>
    !>
    !>  \example{setBracketMin}
    !>  \include{lineno} example/pm_optimization/setBracketMin/main.F90
    !>  \compilef{setBracketMin}
    !>  \output{setBracketMin}
    !>  \include{lineno} example/pm_optimization/setBracketMin/main.out.F90
    !>
    !>  \test
    !>  [test_pm_optimization](@ref test_pm_optimization)
    !>
    !>  \final{setBracketMin}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setBracketMin

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module subroutine setBracketMin_RK5(getFunc, niter, xmin, xlow, xupp, fmin, flow, fupp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBracketMin_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                            :: getFunc
        integer(IK)         , intent(inout)             :: niter
        real(RKG)           , intent(inout)             :: xlow, xupp
        real(RKG)           , intent(out)               :: xmin, fmin
        real(RKG)           , intent(out)   , optional  :: flow, fupp
    end subroutine
#endif

#if RK4_ENABLED
    recursive module subroutine setBracketMin_RK4(getFunc, niter, xmin, xlow, xupp, fmin, flow, fupp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBracketMin_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                            :: getFunc
        integer(IK)         , intent(inout)             :: niter
        real(RKG)           , intent(inout)             :: xlow, xupp
        real(RKG)           , intent(out)               :: xmin, fmin
        real(RKG)           , intent(out)   , optional  :: flow, fupp
    end subroutine
#endif

#if RK3_ENABLED
    recursive module subroutine setBracketMin_RK3(getFunc, niter, xmin, xlow, xupp, fmin, flow, fupp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBracketMin_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                            :: getFunc
        integer(IK)         , intent(inout)             :: niter
        real(RKG)           , intent(inout)             :: xlow, xupp
        real(RKG)           , intent(out)               :: xmin, fmin
        real(RKG)           , intent(out)   , optional  :: flow, fupp
    end subroutine
#endif

#if RK2_ENABLED
    recursive module subroutine setBracketMin_RK2(getFunc, niter, xmin, xlow, xupp, fmin, flow, fupp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBracketMin_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                            :: getFunc
        integer(IK)         , intent(inout)             :: niter
        real(RKG)           , intent(inout)             :: xlow, xupp
        real(RKG)           , intent(out)               :: xmin, fmin
        real(RKG)           , intent(out)   , optional  :: flow, fupp
    end subroutine
#endif

#if RK1_ENABLED
    recursive module subroutine setBracketMin_RK1(getFunc, niter, xmin, xlow, xupp, fmin, flow, fupp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBracketMin_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                            :: getFunc
        integer(IK)         , intent(inout)             :: niter
        real(RKG)           , intent(inout)             :: xlow, xupp
        real(RKG)           , intent(out)               :: xmin, fmin
        real(RKG)           , intent(out)   , optional  :: flow, fupp
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Generate and return the minimum value and the corresponding abscissa `xmin` of the input 1-dimensional function
    !> isolated to a fractional precision of about `tol` using the [Brent method](@ref pm_optimization).<br>
    !>
    !>  \param[in]      getFunc :   The scalar function be minimized.
    !>                              <ol>
    !>                                  <li>    On input, it must take a scalar of the same type and kind as the output argument `xmin`.<br>
    !>                                  <li>    On output, it must return a scalar of the same type and kind as the output argument `xmin`,
    !>                                          containing the function value at the specified input scalar point.<br>
    !>                              </ol>
    !>                              The following demonstrates the interface of `getFunc`,
    !>                              \code{.F90}
    !>                                  function getFunc(x) result(func)
    !>                                      real(RKG), intent(in) :: x
    !>                                      real(RKG) :: func
    !>                                  end function
    !>                              \endcode
    !>                              where `RKG` can refer to any `real` type kind parameter \RKALL supported by the library.<br>
    !>  \param[in]      xlow    :   The input scalar of the same type and kind as the output argument `xmin`,
    !>                              containing the lower bound of the search interval for the function minimum abscissa.<br>
    !>                              The condition `xlow < xmin .and. getFunc(xmin) < getFunc(xlow)` must hold for the corresponding input arguments.<br>
    !>                              (**optional**, default = `0.1`. This requires `0.1` to be within the support of `getFunc()`)
    !>  \param[in]      xupp    :   The input scalar of the same type and kind as the output argument `xmin`,
    !>                              containing the upper bound of the search interval for the function minimum abscissa.<br>
    !>                              The condition `xmin < xupp .and. getFunc(xmin) < getFunc(xupp)` must hold for the corresponding input arguments.<br>
    !>                              (**optional**, default = `0.9`. This requires `0.9` to be within the support of `getFunc()`)
    !>  \param[out]     fmin    :   The output scalar of the same type and kind as the output argument `xmin`.<br>
    !>                              On output, it contains the function value at the identified minimum `xmin`.<br>
    !>                              (**optional**. If missing, the function value at the minimum will not be returned.)
    !>  \param[in]      tol     :   The input positive scalar of the same type and kind as the output argument `xmin`,
    !>                              containing the fractional minimum distance that a new function evaluation point `xmin` can have from any previously evaluated point.<br>
    !>                              Values smaller than the suggestion below might lead to algorithm failure due to roundoff error accumulation.<br>
    !>                              (**optional**. default = `sqrt(epsilon(xmin))`.)
    !>  \param[inout]   niter   :   The input/output positive scalar of type `integer` of default kind \IK containing
    !>                              the **maximum** number of allowed iterations in the algorithm in search of the minimum.<br>
    !>                              On output,
    !>                              <ol>
    !>                                  <li>    **If the algorithm succeeds**, `niter` will be set to the actual number of iterations
    !>                                          taken to find the minimum which is by definition, less than or equal to the input value.<br>
    !>                                  <li>    **If the algorithm fails to converge**, it will be a number larger than the input value (by only `1` unit).<br>
    !>                              </ol>
    !>                              The value of `niter` is effectively the number of calls to the user-specified input function.<br>
    !>                              (**optional**. default = `int(100 * precision(xmin) / 53.)`.)
    !>
    !>  \return
    !>  `xmin`                  :   The output scalar of,
    !>                              <ol>
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the inferred abscissa at the function minimum, if the algorithm succeeds.<br>
    !>
    !>  \interface{getMinBrent}
    !>  \code{.F90}
    !>
    !>      use pm_optimization, only: getMinBrent
    !>
    !>      xmin = getMinBrent(getFunc, xlow = xlow, xupp = xupp, fmin = fmin, tol = tol, niter = niter)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < niter` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < tol` must hold for the corresponding input arguments.<br>
    !>
    !>  \impure
    !>
    !>  \recursive
    !>
    !>  \see
    !>  [getMinBrent](@ref pm_optimization::getMinBrent)<br>
    !>
    !>  \example{getMinBrent}
    !>  \include{lineno} example/pm_optimization/getMinBrent/main.F90
    !>  \compilef{getMinBrent}
    !>  \output{getMinBrent}
    !>  \include{lineno} example/pm_optimization/getMinBrent/main.out.F90
    !>
    !>  \test
    !>  [test_pm_optimization](@ref test_pm_optimization)
    !>
    !>  \final{getMinBrent}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getMinBrent

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module function getMinBrent_RK5(getFunc, xlow, xupp, fmin, tol, niter) result(xmin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinBrent_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                            :: getFunc
        integer(IK)         , intent(inout) , optional  :: niter
        real(RKG)           , intent(in)    , optional  :: xlow, xupp, tol
        real(RKG)           , intent(out)   , optional  :: fmin
        real(RKG)                                       :: xmin
    end function
#endif

#if RK4_ENABLED
    recursive module function getMinBrent_RK4(getFunc, xlow, xupp, fmin, tol, niter) result(xmin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinBrent_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                            :: getFunc
        integer(IK)         , intent(inout) , optional  :: niter
        real(RKG)           , intent(in)    , optional  :: xlow, xupp, tol
        real(RKG)           , intent(out)   , optional  :: fmin
        real(RKG)                                       :: xmin
    end function
#endif

#if RK3_ENABLED
    recursive module function getMinBrent_RK3(getFunc, xlow, xupp, fmin, tol, niter) result(xmin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinBrent_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                            :: getFunc
        integer(IK)         , intent(inout) , optional  :: niter
        real(RKG)           , intent(in)    , optional  :: xlow, xupp, tol
        real(RKG)           , intent(out)   , optional  :: fmin
        real(RKG)                                       :: xmin
    end function
#endif

#if RK2_ENABLED
    recursive module function getMinBrent_RK2(getFunc, xlow, xupp, fmin, tol, niter) result(xmin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinBrent_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                            :: getFunc
        integer(IK)         , intent(inout) , optional  :: niter
        real(RKG)           , intent(in)    , optional  :: xlow, xupp, tol
        real(RKG)           , intent(out)   , optional  :: fmin
        real(RKG)                                       :: xmin
    end function
#endif

#if RK1_ENABLED
    recursive module function getMinBrent_RK1(getFunc, xlow, xupp, fmin, tol, niter) result(xmin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinBrent_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                            :: getFunc
        integer(IK)         , intent(inout) , optional  :: niter
        real(RKG)           , intent(in)    , optional  :: xlow, xupp, tol
        real(RKG)           , intent(out)   , optional  :: fmin
        real(RKG)                                       :: xmin
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Compute and return the minimum value and the corresponding abscissa `xmin` of the input 1-dimensional function
    !> isolated to a fractional precision of about `tol` using the [Brent method](@ref pm_optimization).<br>
    !>
    !>  \param[in]      getFunc :   The scalar function be minimized.
    !>                              <ol>
    !>                                  <li>    On input, it must take a scalar of the same type and kind as the input/output argument `xmin`.<br>
    !>                                  <li>    On output, it must return a scalar of the same type and kind as the input/output argument `xmin`,
    !>                                          containing the function value at the specified input scalar point.<br>
    !>                              </ol>
    !>                              The following demonstrates the interface of `getFunc`,
    !>                              \code{.F90}
    !>                                  function getFunc(x) result(func)
    !>                                      real(RKG), intent(in) :: x
    !>                                      real(RKG) :: func
    !>                                  end function
    !>                              \endcode
    !>                              where `RKG` can refer to any `real` type kind parameter \RKALL supported by the library.<br>
    !>  \param[inout]   xmin    :   The input/output scalar of type `real` of kind \RKALL.<br>
    !>                              <ol>
    !>                                  <li>    On input, it must contain the initial best guess abscissa at the function minimum.<br>
    !>                                          Note that the specified input value **must be between** the input arguments `xlow` and `xupp`,
    !>                                          such that the triplet `(xlow, getFunc(xlow)), (xmin, getFunc(fmin)), (xupp, getFunc(fupp))` can be interpolated by a convex parabolic curve.<br>
    !>                                  <li>    On output, it will contain the inferred abscissa at the function minimum, if the algorithm succeeds.<br>
    !>                              </ol>
    !>  \param[in]      xlow    :   The input scalar of the same type and kind as the input/output argument `xmin`,
    !>                              containing the lower bound of the search interval for the function minimum abscissa.<br>
    !>                              The condition `xlow < xmin .and. getFunc(xmin) < getFunc(xlow)` must hold for the corresponding input arguments.<br>
    !>  \param[in]      xupp    :   The input scalar of the same type and kind as the input/output argument `xmin`,
    !>                              containing the upper bound of the search interval for the function minimum abscissa.<br>
    !>                              The condition `xmin < xupp .and. getFunc(xmin) < getFunc(xupp)` must hold for the corresponding input arguments.<br>
    !>  \param[inout]   fmin    :   The input/output scalar of the same type and kind as the input/output argument `xmin`.<br>
    !>                              <ol>
    !>                                  <li>    On input, it must contain `getFunc(xmin)`.<br>
    !>                                  <li>    On output, it will contain the function value at the identified minimum `xmin`.<br>
    !>                              </ol>
    !>  \param[in]      tol     :   The input positive scalar of the same type and kind as the input/output argument `xmin`,
    !>                              containing the minimum distance that a new function evaluation point `xmin` can have from any previously evaluated point.<br>
    !>                              Values smaller than the suggestion below might lead to algorithm failure due to roundoff error accumulation.<br>
    !>                              A reasonable choice is `tol = sqrt(epsilon(1._RKG))`.<br>
    !>  \param[inout]   niter   :   The input/output positive scalar of type `integer` of default kind \IK containing
    !>                              the **maximum** number of allowed iterations in the algorithm in search of the minimum.<br>
    !>                              On output,
    !>                              <ol>
    !>                                  <li>    **If the algorithm succeeds**, `niter` will be set to the actual number of iterations
    !>                                          taken to find the minimum which is by definition, less than or equal to the input value.<br>
    !>                                  <li>    **If the algorithm fails to converge**, it will be a number larger than the input value (by only `1` unit).<br>
    !>                              </ol>
    !>                              The value of `niter` is effectively the number of calls to the user-specified input function.<br>
    !>                              A reasonable choice is `niter = int(100 * precision(xmin) / 53.)`.<br>
    !>
    !>  \interface{setMinBrent}
    !>  \code{.F90}
    !>
    !>      use pm_optimization, only: setMinBrent
    !>
    !>      call setMinBrent(getFunc, xmin, xlow, xupp, fmin, tol, niter)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `xlow < xmin .and. xmin < xupp` must hold for the corresponding input arguments.<br>
    !>  The condition `getFunc(xlow) > fmin .or. fmin < getFunc(xupp)` must hold for the corresponding input arguments.<br>
    !>  The condition `fmin == getFunc(xmin)` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < niter` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < tol` must hold for the corresponding input arguments.<br>
    !>
    !>  \impure
    !>
    !>  \recursive
    !>
    !>  \see
    !>  [getMinBrent](@ref pm_optimization::getMinBrent)<br>
    !>  [setMinBrent](@ref pm_optimization::setMinBrent)<br>
    !>
    !>  \example{setMinBrent}
    !>  \include{lineno} example/pm_optimization/setMinBrent/main.F90
    !>  \compilef{setMinBrent}
    !>  \output{setMinBrent}
    !>  \include{lineno} example/pm_optimization/setMinBrent/main.out.F90
    !>
    !>  \test
    !>  [test_pm_optimization](@ref test_pm_optimization)
    !>
    !>  \final{setMinBrent}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setMinBrent

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module subroutine setMinBrent_RK5(getFunc, xmin, xlow, xupp, fmin, tol, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinBrent_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                            :: getFunc
        real(RKG)           , value                     :: xlow, xupp, tol
        real(RKG)           , intent(inout)             :: fmin, xmin
        integer(IK)         , intent(inout)             :: niter
    end subroutine
#endif

#if RK4_ENABLED
    recursive module subroutine setMinBrent_RK4(getFunc, xmin, xlow, xupp, fmin, tol, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinBrent_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                            :: getFunc
        real(RKG)           , value                     :: xlow, xupp, tol
        real(RKG)           , intent(inout)             :: fmin, xmin
        integer(IK)         , intent(inout)             :: niter
    end subroutine
#endif

#if RK3_ENABLED
    recursive module subroutine setMinBrent_RK3(getFunc, xmin, xlow, xupp, fmin, tol, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinBrent_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                            :: getFunc
        real(RKG)           , value                     :: xlow, xupp, tol
        real(RKG)           , intent(inout)             :: fmin, xmin
        integer(IK)         , intent(inout)             :: niter
    end subroutine
#endif

#if RK2_ENABLED
    recursive module subroutine setMinBrent_RK2(getFunc, xmin, xlow, xupp, fmin, tol, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinBrent_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                            :: getFunc
        real(RKG)           , value                     :: xlow, xupp, tol
        real(RKG)           , intent(inout)             :: fmin, xmin
        integer(IK)         , intent(inout)             :: niter
    end subroutine
#endif

#if RK1_ENABLED
    recursive module subroutine setMinBrent_RK1(getFunc, xmin, xlow, xupp, fmin, tol, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinBrent_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                            :: getFunc
        real(RKG)           , value                     :: xlow, xupp, tol
        real(RKG)           , intent(inout)             :: fmin, xmin
        integer(IK)         , intent(inout)             :: niter
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Generate and return `.true.` **if and only if** the algorithm **fails** to find the minimum value and the corresponding abscissa `xmin(1:ndim)` of the input arbitrary (`ndim`) dimensional-support
    !> function isolated to a fractional precision of about `tol` using the [Powell unconstrained derivative-free minimization method](@ref pm_optimization).<br>
    !>
    !>  \param[in]      getFunc :   The scalar function be minimized.
    !>                              <ol>
    !>                                  <li>    On input, it must take a vector of size `ndim` of the same type and kind as the input/output argument `xmin`.<br>
    !>                                  <li>    On output, it must return a scalar of the same type and kind as the input/output argument `xmin`,
    !>                                          containing the function value at the specified input vector point.<br>
    !>                              </ol>
    !>                              The following demonstrates the interface of `getFunc`,
    !>                              \code{.F90}
    !>                                  function getFunc(x) result(func)
    !>                                      real(RKG), intent(in) :: x(ndim)
    !>                                      real(RKG) :: func
    !>                                  end function
    !>                              \endcode
    !>                              where `RKG` can refer to any `real` type kind parameter \RKALL supported by the library.<br>
    !>  \param[inout]   xmin    :   The input/output vector of size `ndim` of type `real` of kind \RKALL.<br>
    !>                              <ol>
    !>                                  <li>    On input, it must contain the initial best guess abscissa at the function minimum.<br>
    !>                                  <li>    On output, it will contain the inferred abscissa at the function minimum, if the algorithm succeeds.<br>
    !>                              </ol>
    !>  \param[inout]   fmin    :   The input/output scalar of the same type and kind as the input/output argument `xmin`.<br>
    !>                              <ol>
    !>                                  <li>    On input, it must contain `getFunc(xmin)`.<br>
    !>                                  <li>    On output, it will contain the function value at the identified minimum abscissa `xmin`.<br>
    !>                              </ol>
    !>                              (**optional**, default = `getFunc(xmin)`)
    !>  \param[in]      tol     :   The input positive scalar of the same type and kind as the input/output argument `xmin`,
    !>                              containing the minimum distance that a new function evaluation point `xmin` can have from any previously evaluated point.<br>
    !>                              Values smaller than the suggestion below might lead to algorithm failure due to roundoff error accumulation.<br>
    !>                              (**optional**, default = `sqrt(epsilon(1._RKG))`).<br>
    !>  \param[in]      niter   :   The input positive scalar of type `integer` of default kind \IK containing
    !>                              the **maximum** number of allowed iterations at every step of the algorithm in search of a univariate minimum along a specific direction.<br>
    !>                              (**optional**, default = `int(100 * precision(xmin) / 53.)`.)
    !>
    !>  \return
    !>  `failed`                :   The output scalar of type `logical` of default kind \LK that is `.true.` **if and only if** the algorithm
    !>                              fails to find the minimum of the function, otherwise it is `.false.` indicating convergence to the minimum.<br>
    !>
    !>  \interface{isFailedMinPowell}
    !>  \code{.F90}
    !>
    !>      use pm_optimization, only: isFailedMinPowell
    !>      use pm_kind, only: LK
    !>      logical(LK) :: failed
    !>
    !>      failed = isFailedMinPowell(getFunc, xmin(1:ndim), fmin = fmin, tol = tol, niter = niter)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `fmin == getFunc(xmin)` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < niter` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < tol` must hold for the corresponding input arguments.<br>
    !>
    !>  \impure
    !>
    !>  \recursive
    !>
    !>  \see
    !>  [getMinBrent](@ref pm_optimization::getMinBrent)<br>
    !>  [setMinBrent](@ref pm_optimization::setMinBrent)<br>
    !>  [isFailedMinPowell](@ref pm_optimization::isFailedMinPowell)<br>
    !>  [setMinPowell](@ref pm_optimization::setMinPowell)<br>
    !>
    !>  \example{isFailedMinPowell}
    !>  \include{lineno} example/pm_optimization/isFailedMinPowell/main.F90
    !>  \compilef{isFailedMinPowell}
    !>  \output{isFailedMinPowell}
    !>  \include{lineno} example/pm_optimization/isFailedMinPowell/main.out.F90
    !>
    !>  \test
    !>  [test_pm_optimization](@ref test_pm_optimization)
    !>
    !>  \final{isFailedMinPowell}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isFailedMinPowell

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module function isFailedMinPowell_RK5(getFunc, xmin, fmin, tol, niter) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedMinPowell_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                            :: getFunc
        real(RKG)           , intent(inout) , contiguous                :: xmin(:)
        real(RKG)           , intent(inout) , optional                  :: fmin
        real(RKG)           , intent(in)    , optional                  :: tol
        integer(IK)         , intent(in)    , optional                  :: niter
        logical(LK)                                                     :: failed
    end function
#endif

#if RK4_ENABLED
    recursive module function isFailedMinPowell_RK4(getFunc, xmin, fmin, tol, niter) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedMinPowell_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                            :: getFunc
        real(RKG)           , intent(inout) , contiguous                :: xmin(:)
        real(RKG)           , intent(inout) , optional                  :: fmin
        real(RKG)           , intent(in)    , optional                  :: tol
        integer(IK)         , intent(in)    , optional                  :: niter
        logical(LK)                                                     :: failed
    end function
#endif

#if RK3_ENABLED
    recursive module function isFailedMinPowell_RK3(getFunc, xmin, fmin, tol, niter) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedMinPowell_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                            :: getFunc
        real(RKG)           , intent(inout) , contiguous                :: xmin(:)
        real(RKG)           , intent(inout) , optional                  :: fmin
        real(RKG)           , intent(in)    , optional                  :: tol
        integer(IK)         , intent(in)    , optional                  :: niter
        logical(LK)                                                     :: failed
    end function
#endif

#if RK2_ENABLED
    recursive module function isFailedMinPowell_RK2(getFunc, xmin, fmin, tol, niter) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedMinPowell_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                            :: getFunc
        real(RKG)           , intent(inout) , contiguous                :: xmin(:)
        real(RKG)           , intent(inout) , optional                  :: fmin
        real(RKG)           , intent(in)    , optional                  :: tol
        integer(IK)         , intent(in)    , optional                  :: niter
        logical(LK)                                                     :: failed
    end function
#endif

#if RK1_ENABLED
    recursive module function isFailedMinPowell_RK1(getFunc, xmin, fmin, tol, niter) result(failed)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isFailedMinPowell_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                            :: getFunc
        real(RKG)           , intent(inout) , contiguous                :: xmin(:)
        real(RKG)           , intent(inout) , optional                  :: fmin
        real(RKG)           , intent(in)    , optional                  :: tol
        integer(IK)         , intent(in)    , optional                  :: niter
        logical(LK)                                                     :: failed
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Compute and return the minimum value and the corresponding abscissa `xmin(1:ndim)` of the input arbitrary (`ndim`) dimensional-support
    !> function isolated to a fractional precision of about `tol` using the [Powell unconstrained derivative-free minimization method](@ref pm_optimization).<br>
    !>
    !>  \param[in]      getFunc :   The scalar function be minimized.
    !>                              <ol>
    !>                                  <li>    On input, it must take a vector of size `ndim` of the same type and kind as the input/output argument `xmin`.<br>
    !>                                  <li>    On output, it must return a scalar of the same type and kind as the input/output argument `xmin`,
    !>                                          containing the function value at the specified input vector point.<br>
    !>                              </ol>
    !>                              The following demonstrates the interface of `getFunc`,
    !>                              \code{.F90}
    !>                                  function getFunc(x) result(func)
    !>                                      real(RKG), intent(in) :: x(ndim)
    !>                                      real(RKG) :: func
    !>                                  end function
    !>                              \endcode
    !>                              where `RKG` can refer to any `real` type kind parameter \RKALL supported by the library.<br>
    !>  \param[inout]   xmin    :   The input/output vector of size `ndim` of type `real` of kind \RKALL.<br>
    !>                              <ol>
    !>                                  <li>    On input, it must contain the initial best guess abscissa at the function minimum.<br>
    !>                                  <li>    On output, it will contain the inferred abscissa at the function minimum, if the algorithm succeeds.<br>
    !>                              </ol>
    !>  \param[inout]   fmin    :   The input/output scalar of the same type and kind as the input/output argument `xmin`.<br>
    !>                              <ol>
    !>                                  <li>    On input, it must contain `getFunc(xmin)`.<br>
    !>                                  <li>    On output, it will contain the function value at the identified minimum abscissa `xmin`.<br>
    !>                              </ol>
    !>  \param[inout]   dirset  :   The input/output matrix of the same type and kind as the input/output argument `xmin`,
    !>                              <ol>
    !>                                  <li>    On input, it must contain the initial **direction set** matrix, each column of which is a unit vector along the corresponding coordinates axis.<br>
    !>                                  <li>    On output, the initial contents will be destroyed because the matrix is used as a workspace.<br>
    !>                              </ol>
    !>  \param[in]      tol     :   The input positive scalar of the same type and kind as the input/output argument `xmin`,
    !>                              containing the minimum distance that a new function evaluation point `xmin` can have from any previously evaluated point.<br>
    !>                              Values smaller than the suggestion below might lead to algorithm failure due to roundoff error accumulation.<br>
    !>                              A reasonable choice is `tol = sqrt(epsilon(1._RKG))`.<br>
    !>  \param[inout]   niter   :   The input/output positive scalar of type `integer` of default kind \IK containing
    !>                              the **maximum** number of allowed iterations at every step of the algorithm in search of a univariate minimum along a specific direction.<br>
    !>                              On output,
    !>                              <ol>
    !>                                  <li>    **If the algorithm succeeds**, `niter` will be set to the actual number of iterations
    !>                                          taken to find the minimum which is by definition, less than or equal to the input value.<br>
    !>                                  <li>    **If the algorithm fails to converge**, it will be a number larger than the input value (by only `1` unit).<br>
    !>                              </ol>
    !>                              A reasonable choice is `niter = int(100 * precision(xmin) / 53.)`.<br>
    !>
    !>  \interface{setMinPowell}
    !>  \code{.F90}
    !>
    !>      use pm_optimization, only: setMinPowell
    !>
    !>      call setMinPowell(getFunc, xmin, fmin, dirset, tol, niter)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `fmin == getFunc(xmin)` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < niter` must hold for the corresponding input arguments.<br>
    !>  The condition `0 < tol` must hold for the corresponding input arguments.<br>
    !>
    !>  \impure
    !>
    !>  \recursive
    !>
    !>  \see
    !>  [getMinBrent](@ref pm_optimization::getMinBrent)<br>
    !>  [setMinBrent](@ref pm_optimization::setMinBrent)<br>
    !>  [isFailedMinPowell](@ref pm_optimization::isFailedMinPowell)<br>
    !>  [setMinPowell](@ref pm_optimization::setMinPowell)<br>
    !>
    !>  \example{setMinPowell}
    !>  \include{lineno} example/pm_optimization/setMinPowell/main.F90
    !>  \compilef{setMinPowell}
    !>  \output{setMinPowell}
    !>  \include{lineno} example/pm_optimization/setMinPowell/main.out.F90
    !>
    !>  \test
    !>  [test_pm_optimization](@ref test_pm_optimization)
    !>
    !>  \final{setMinPowell}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday March 7, 2017, 3:50 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setMinPowell

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module subroutine setMinPowell_RK5(getFunc, xmin, fmin, dirset, tol, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinPowell_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                                :: getFunc
        real(RKG)           , intent(inout) , contiguous    :: xmin(:), dirset(:,:)
        real(RKG)           , intent(inout)                 :: fmin
        real(RKG)           , intent(in)                    :: tol
        integer(IK)         , intent(inout)                 :: niter
    end subroutine
#endif

#if RK4_ENABLED
    recursive module subroutine setMinPowell_RK4(getFunc, xmin, fmin, dirset, tol, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinPowell_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                                :: getFunc
        real(RKG)           , intent(inout) , contiguous    :: xmin(:), dirset(:,:)
        real(RKG)           , intent(inout)                 :: fmin
        real(RKG)           , intent(in)                    :: tol
        integer(IK)         , intent(inout)                 :: niter
    end subroutine
#endif

#if RK3_ENABLED
    recursive module subroutine setMinPowell_RK3(getFunc, xmin, fmin, dirset, tol, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinPowell_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                                :: getFunc
        real(RKG)           , intent(inout) , contiguous    :: xmin(:), dirset(:,:)
        real(RKG)           , intent(inout)                 :: fmin
        real(RKG)           , intent(in)                    :: tol
        integer(IK)         , intent(inout)                 :: niter
    end subroutine
#endif

#if RK2_ENABLED
    recursive module subroutine setMinPowell_RK2(getFunc, xmin, fmin, dirset, tol, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinPowell_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                                :: getFunc
        real(RKG)           , intent(inout) , contiguous    :: xmin(:), dirset(:,:)
        real(RKG)           , intent(inout)                 :: fmin
        real(RKG)           , intent(in)                    :: tol
        integer(IK)         , intent(inout)                 :: niter
    end subroutine
#endif

#if RK1_ENABLED
    recursive module subroutine setMinPowell_RK1(getFunc, xmin, fmin, dirset, tol, niter)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMinPowell_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                                :: getFunc
        real(RKG)           , intent(inout) , contiguous    :: xmin(:), dirset(:,:)
        real(RKG)           , intent(inout)                 :: fmin
        real(RKG)           , intent(in)                    :: tol
        integer(IK)         , intent(inout)                 :: niter
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_optimization ! LCOV_EXCL_LINE