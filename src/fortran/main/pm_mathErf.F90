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
!>  This module contains classes and procedures for computing the mathematical Inverse Error Function.
!>
!>  \details
!>  The **error function** (or the **Gauss error function**),
!>  denoted by \f$\ms{erf}(\cdot)\f$, is a complex function of a complex variable defined as,
!>  \f{equation}{
!>      \ms{erf}(x) = {\frac{2}{\sqrt{\pi}}} \int_{0}^{x} e^{-t^{2}} ~ \mathrm{d}t ~.
!>  \f}
!>  This integral is a special (non-elementary) Sigmoid function that
!>  occurs in probability, statistics, and partial differential equations.<br>
!>  In many of these applications, the function argument is however a real number.<br>
!>  If the function argument is real, then the function value is also real.<br>
!>  For non-negative values of \f$x\f$, the error function has the following interpretation:
!>  For a random variable \f$Y\f$ that is [normally distributed](@ref pm_distNorm)
!>  with mean \f$0\f$ and standard deviation \f$\frac{1}{\sqrt{2}}\f$,
!>  \f$\ms{erf}(x)\f$ is the probability that \f$Y\f$ falls in the range \f$[−x, x]\f$.<br>
!>  Two closely related functions are the **complementary error function (erfc)** defined as,
!>  \f{equation}{
!>      \ms{erfc}(x) = 1 - \ms{erf}(x) ~,
!>  \f}
!>  and the **imaginary error function (erfi)** defined as,
!>  \f{equation}{
!>      \ms{erfi}(x) = -i\ms{erf}(ix) ~,
!>  \f}
!>  where \f$i\f$ is the imaginary unit.<br>
!>
!>  <b>Inverse error function</b><br>
!>  Given a complex number \f$x\f$, there is not a unique complex number \f$w\f$ satisfying \f$\ms{erf}(w) = x\f$.<br>
!>  Therefore, a true inverse function would be multivalued.<br>
!>  However, for \f$−1 < x < 1\f$, there is a unique real number denoted \f$\ms{erf}^{−1}(x)\f$ satisfying,
!>  \f{equation}{
!>      \ms{erf}\left(\ms{erf}^{-1}(x)\right) = x ~.
!>  \f}
!>  The inverse error function is usually defined with domain \f$(−1, 1)\f$
!>  and it is restricted to this domain in many computer algebra systems.<br>
!>  However, it can be extended to the disk \f$|x| < 1\f$ of the complex plane, using the Maclaurin series,
!>  \f{equation}{
!>      \ms{erf}^{-1}(x) = \sum_{k=0}^{\infty} \frac{c_{k}}{2k+1} \left({\frac{\sqrt{\pi}}{2}} x \right)^{2k+1} ~,
!>  \f}
!>  where \f$c_0 = 1\f$ and,
!>  \f{equation}{
!>      \begin{aligned}
!>          c_k & = \sum_{m=0}^{k-1} \frac{c_{m} c_{k-1-m}}{(m+1) (2m+1)} \\
!>              & = \left\{1, 1, \frac{7}{6}, \frac{127}{90}, \frac{4369}{2520}, \frac{34807}{16200}, \ldots \right\} ~.
!>      \end{aligned}
!>  \f}
!>  Therefore,
!>  \f{equation}{
!>      \ms{erf}^{-1}(x) = \frac{\sqrt{\pi}}{2} \left(x + \frac{\pi}{12} x^{3} + \frac{7\pi^{2}}{480} x^{5} + \frac{127\pi^{3}}{40320} x^{7} + \frac{4369\pi^{4}}{5806080} x^{9} + \frac{34807\pi^{5}}{182476800} x^{11} + \cdots \right) ~.
!>  \f}
!>  The error function value at \f$x = \pm\infty\f$ is equal to \f$\pm1\f$.<br>
!>  For \f$|x| < 1\f$, \f$\ms{erf}\left(\ms{erf}^{−1}(x)\right) = x\f$.<br>
!>  The **inverse complementary error function** is defined as,
!>  \f{equation}{
!>      \ms{erfc}^{-1}(1 - x) = \ms{erf}^{-1}(x) ~.
!>  \f}
!>  For real \f$x\f$, there is a unique real number \f$\ms{erfi}^{−1}(x)\f$ satisfying \f$\ms{erfi}(\ms{erfi}^{−1}(x) = x\f$.<br>
!>  The **inverse imaginary error function** is defined as \f$\ms{erfi}^{−1}(x)\f$.<br>
!>
!>  <b>Numerical computation</b><br>
!>
!>  The `erf()` and `erfc()` intrinsic Fortran functions readily return the value of
!>  the Error function at any given input value with arbitrary `real` type and kind.<br>
!>  Theoretically, for any real \f$x\f$, the [Newton root-finding method](@ref pm_mathRoot)
!>  can be used to compute the **inverse error function** \f$\ms{erfi}^{−1}(x)\f$
!>  and for −1 ≤ x ≤ 1, the following Maclaurin series converges:
!>  \f{equation}{
!>      \ms{erfi}^{-1}(x) = \sum_{k=0}^{\infty}{\frac{(-1)^{k}c_{k}}{2k+1}}\left({\frac{\sqrt{\pi}}{2}}(x)\right)^{2k+1} ~,
!>  \f}
!>  where \f$c_k\f$ is defined as above.
!>
!>  The procedures of this module combine multiple varying-precision approaches to make a
!>  decision at compile-time about the best strategy for computing the inverse error function.<br>
!>
!>  \see
!>  [pm_mathBeta](@ref pm_mathBeta)<br>
!>  [pm_mathErf](@ref pm_mathErf)<br>
!>  [pm_mathGamma](@ref pm_mathGamma)<br>
!>  [pm_distNorm](@ref pm_distNorm)<br>
!>  [getNormQuan](@ref pm_distNorm::getNormQuan)<br>
!>  [setNormQuan](@ref pm_distNorm::setNormQuan)<br>
!>
!>  \test
!>  [test_pm_mathErf](@ref test_pm_mathErf)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Nov 10, 2009, 8:53 PM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_mathErf

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathErf"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Inverse Error Function \f$\ms{erf}^{-1}(x)\f$ for an input `real` value in range \f$(-1, +1)\f$
    !>  as defined in the details section of [pm_mathErf](@ref pm_mathErf).
    !>
    !>  \details
    !>  The **complementary inverse error function** \f$\ms{erfc}^{-1}(y), y\in(0, 2)\f$ can be readily computed as `getErfInv(1 - y)`.<br>
    !>  See the documentation of [pm_mathErf](@ref pm_mathErf) for more information.<br>
    !>
    !>  \param[in]  x       :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as the output `erfinv`.<br>
    !>  \param[in]  abserr  :   The input positive scalar (or array of the same shape as other array-like arguments) of the same type and kind as the output `erfinv`.<br>
    !>                          If present, the output approximate value for the inverse error function is guaranteed to be at most `abserr` away from the true value \f$\ms{erf}^{-1}(x)\f$.<br>
    !>                          This guarantee is currently tested and validated up to \f$2\times10^{-26}\f$.<br>
    !>                          (**optional**, default = `100 * epsilon(x)`)
    !>  \return
    !>  `erfinv`            :   The output scalar (or array of the same shape as input array-like arguments) of,
    !>                          <ol>
    !>                              <li>    type `real` of kind \RKALL,<br>
    !>                          </ol>
    !>                          containing the approximate value of the inverse error function at the specified input point `x`.
    !>
    !>  \interface{getErfInv}
    !>  \code{.F90}
    !>
    !>      use pm_mathErf, only: getErfInv
    !>
    !>      erfinv = getErfInv(x, abserr = abserr)
    !>      erfinv(..) = getErfInv(x(..), abserr = abserr(..))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < abserr` must hold for the corresponding input arguments.<br>
    !>  The conditions `-1 < x .and. x < 1` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getErfInv](@ref pm_mathErf::getErfInv)<br>
    !>  [setErfInv](@ref pm_mathErf::setErfInv)<br>
    !>  [getNormCDF](@ref pm_distNorm::getNormCDF)<br>
    !>  [setNormCDF](@ref pm_distNorm::setNormCDF)<br>
    !>  [getNormQuan](@ref pm_distNorm::getNormQuan)<br>
    !>  [setNormQuan](@ref pm_distNorm::setNormQuan)<br>
    !>  [getNormRand](@ref pm_distNorm::getNormRand)<br>
    !>  [setNormRand](@ref pm_distNorm::setNormRand)<br>
    !>  [getNormLogPDF](@ref pm_distNorm::getNormLogPDF)<br>
    !>  [setNormLogPDF](@ref pm_distNorm::setNormLogPDF)<br>
    !>
    !>  \example{getErfInv}
    !>  \include{lineno} example/pm_mathErf/getErfInv/main.F90
    !>  \compilef{getErfInv}
    !>  \output{getErfInv}
    !>  \include{lineno} example/pm_mathErf/getErfInv/main.out.F90
    !>  \postproc{getErfInv}
    !>  \include{lineno} example/pm_mathErf/getErfInv/main.py
    !>  \vis{getErfInv}
    !>  \image html pm_mathErf/getErfInv/getErfInv.RK.png width=700
    !>  \image html pm_mathErf/getErfInv/getErfInv.RKS.abserr.png width=700
    !>  \image html pm_mathErf/getErfInv/getErfInv.RKD.abserr.png width=700
    !>  \image html pm_mathErf/getErfInv/getErfInv.RKH.abserr.png width=700
    !>
    !>  \test
    !>  [test_pm_mathErf](@ref test_pm_mathErf)
    !>
    !>  \final{getErfInv}
    !>
    !>  \author
    !>  \AmirShahmoradi, Nov 10, 2009, 8:53 PM, Michigan
    interface getErfInv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getErfInv_RK5(x, abserr) result(erfinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getErfInv_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)            :: x
        real(RKG)   , intent(in), optional  :: abserr
        real(RKG)                           :: erfinv
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getErfInv_RK4(x, abserr) result(erfinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getErfInv_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)            :: x
        real(RKG)   , intent(in), optional  :: abserr
        real(RKG)                           :: erfinv
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getErfInv_RK3(x, abserr) result(erfinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getErfInv_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)            :: x
        real(RKG)   , intent(in), optional  :: abserr
        real(RKG)                           :: erfinv
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getErfInv_RK2(x, abserr) result(erfinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getErfInv_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)            :: x
        real(RKG)   , intent(in), optional  :: abserr
        real(RKG)                           :: erfinv
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getErfInv_RK1(x, abserr) result(erfinv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getErfInv_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)            :: x
        real(RKG)   , intent(in), optional  :: abserr
        real(RKG)                           :: erfinv
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Inverse Error Function \f$\ms{erf}^{-1}(x)\f$ for an input `real` value in range \f$(-1, +1)\f$
    !>  as defined in the details section of [pm_mathErf](@ref pm_mathErf).
    !>
    !>  \details
    !>  The **complementary inverse error function** \f$\ms{erfc}^{-1}(y), y\in(0, 2)\f$ can be readily computed as `call setErfInv(1 - y)`.<br>
    !>  See the documentation of [pm_mathErf](@ref pm_mathErf) for more information.<br>
    !>
    !>  \param[in]  erfinv  :   The output scalar (or array of the same shape as input array-like arguments) of,
    !>                          <ol>
    !>                              <li>    type `real` of kind \RKALL,<br>
    !>                          </ol>
    !>                          containing the approximate value of the inverse error function at the specified input point `x`.
    !>  \param[in]  x       :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as the output `erfinv`.<br>
    !>  \param[in]  abserr  :   The input positive scalar (or array of the same shape as other array-like arguments) of the same type and kind as the output `erfinv`.<br>
    !>                          If present, the output approximate value for the inverse error function is guaranteed to be at most `abserr` away from the true value \f$\ms{erf}^{-1}(x)\f$.<br>
    !>                          This guarantee is currently tested and validated up to \f$2\times10^{-26}\f$.<br>
    !>                          A reasonable value could be `epsilon(x)**0.66`.<br>
    !>                          The current implementation of this generic interface contains three methods of
    !>                          computing the inverse error function corresponding to three levels of increasing accuracy,
    !>                          <ol>
    !>                              <li>    For any `1.e-7 < abserr`, the output `erfinv` will have at most an absolute error `1.e-7` with respect to theoretical value.<br>
    !>                              <li>    For any `2.e-24 < abserr`, the output `erfinv` will have at most an absolute error `2.e-24` with respect to theoretical value.<br>
    !>                              <li>    For any `abserr < 2.e-24`, the output `erfinv` is currently verified to have at most an absolute error `2.e-26` with respect to theoretical value.<br>
    !>                                      For a wide range of input `x` values, the absolute error is practically orders of magnitude smaller than this verified upper limit for the error.<br>
    !>                                      See the example illustrations below for more information.<br>
    !>                          </ol>
    !>
    !>  \interface{setErfInv}
    !>  \code{.F90}
    !>
    !>      use pm_mathErf, only: setErfInv
    !>
    !>      call setErfInv(erfinv, x, abserr)
    !>      call setErfInv(erfinv, x(..), abserr(..))
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < abserr` must hold for the corresponding input arguments.<br>
    !>  The conditions `-1 < x .and. x < 1` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getErfInv](@ref pm_mathErf::getErfInv)<br>
    !>  [setErfInv](@ref pm_mathErf::setErfInv)<br>
    !>  [getNormCDF](@ref pm_distNorm::getNormCDF)<br>
    !>  [setNormCDF](@ref pm_distNorm::setNormCDF)<br>
    !>  [getNormQuan](@ref pm_distNorm::getNormQuan)<br>
    !>  [setNormQuan](@ref pm_distNorm::setNormQuan)<br>
    !>  [getNormRand](@ref pm_distNorm::getNormRand)<br>
    !>  [setNormRand](@ref pm_distNorm::setNormRand)<br>
    !>  [getNormLogPDF](@ref pm_distNorm::getNormLogPDF)<br>
    !>  [setNormLogPDF](@ref pm_distNorm::setNormLogPDF)<br>
    !>
    !>  \example{setErfInv}
    !>  \include{lineno} example/pm_mathErf/setErfInv/main.F90
    !>  \compilef{setErfInv}
    !>  \output{setErfInv}
    !>  \include{lineno} example/pm_mathErf/setErfInv/main.out.F90
    !>  \postproc{setErfInv}
    !>  \include{lineno} example/pm_mathErf/setErfInv/main.py
    !>  \vis{setErfInv}
    !>  \image html pm_mathErf/setErfInv/setErfInv.RK.png width=700
    !>  \image html pm_mathErf/setErfInv/setErfInv.RKS.abserr.png width=700
    !>  \image html pm_mathErf/setErfInv/setErfInv.RKD.abserr.png width=700
    !>  \image html pm_mathErf/setErfInv/setErfInv.RKH.abserr.png width=700
    !>
    !>  \test
    !>  [test_pm_mathErf](@ref test_pm_mathErf)
    !>
    !>  \final{setErfInv}
    !>
    !>  \author
    !>  \AmirShahmoradi, Nov 10, 2009, 8:53 PM, Michigan
    interface setErfInv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setErfInv_RK5(erfinv, x, abserr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setErfInv_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)   :: erfinv
        real(RKG)   , intent(in)    :: x, abserr
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setErfInv_RK4(erfinv, x, abserr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setErfInv_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)   :: erfinv
        real(RKG)   , intent(in)    :: x, abserr
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setErfInv_RK3(erfinv, x, abserr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setErfInv_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)   :: erfinv
        real(RKG)   , intent(in)    :: x, abserr
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setErfInv_RK2(erfinv, x, abserr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setErfInv_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)   :: erfinv
        real(RKG)   , intent(in)    :: x, abserr
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setErfInv_RK1(erfinv, x, abserr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setErfInv_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)   :: erfinv
        real(RKG)   , intent(in)    :: x, abserr
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_mathErf