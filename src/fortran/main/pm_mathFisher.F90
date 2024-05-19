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
!>  This module contains procedures and generic interfaces for evaluating the Fisher transformation and its inverse.<br>
!>
!>  \details
!>  The Fisher transformation (or **Fisher z-transformation**) of a Pearson correlation coefficient is its inverse hyperbolic tangent (`atanh`).<br>
!>  When the **sample correlation coefficient** \f$r\f$ is near \f$+1\f$ or \f$-1\f$, its distribution is highly skewed, which makes it difficult
!>  to estimate confidence intervals and apply tests of significance for the **population correlation coefficient** \f$\rho\f$.<br>
!>  The Fisher transformation solves this problem by yielding a variable whose distribution is approximately normally distributed,
!>  with a variance that is stable over different values of \f$r\f$.<br>
!>
!>  Definition
!>  ----------
!>
!>  Given a set of \f$N\f$ bivariate sample pairs \f$(X_i, Y_i), i = 1, ..., N\f$, the sample correlation coefficient \f$r\f$ is given by,
!>  \f{equation}{
!>      r = {\frac{\up{cov}(X, Y)}{\sigma_{X}\sigma_{Y}}} =
!>      \frac{
!>          \sum_{i = 1}^{N}(X_{i} - {\bar{X}})(Y_{i} - {\bar {Y}})
!>      }{
!>          {\sqrt{\sum_{i = 1}^{N}(X_{i} - {\bar{X}})^{2}}}{\sqrt{\sum_{i=1}^{N}(Y_{i} - {\bar {Y}})^{2}}}
!>      }
!>  \f}
!>
!>  Here \f$\up{cov}(X, Y)\f$ stands for the covariance between the variables \f$X\f$ and \f$Y\f$ and \f$\sigma\f$ stands
!>  for the standard deviation of the respective variable. The Fisher z-transformation of \f$r\f$ is defined as,
!>  \f{equation}{
!>      z = {1 \over 2} \ln \left( {1 + r \over 1 - r} \right) = \up{atanh}(r) ~,
!>  \f}
!>  where \f$\ln\f$ is the natural logarithm function and `atanh` is the inverse hyperbolic tangent function.<br>
!>
!>  If \f$(X, Y)\f$ has a bivariate normal distribution with correlation \f$\rho\f$ and the pairs \f$(X_i, Y_i)\f$ are independent and identically distributed,
!>  then \f$z\f$ is approximately normally distributed with mean,
!>  \f{equation}{
!>      {1 \over 2} \ln \left({{1 + \rho} \over {1 - \rho}} \right) ~,
!>  \f}
!>
!>  and standard deviation
!>  \f{equation}{
!>      1 \over {\sqrt {N - 3}} ~,
!>  \f}
!>
!>  where \f$N\f$ is the sample size, and \f$\rho\f$ is the true correlation coefficient.<br>
!>
!>  This transformation and its inverse,
!>  \f{equation}{
!>      r = {\frac{\exp(2z) - 1}{\exp(2z) + 1}} = \up{tanh}(z) ~,
!>  \f}
!>  can be used to construct a large-sample confidence interval for \f$r\f$ using standard normal theory and derivations.<br>
!>
!>  \note
!>  The Fisher transformation can be applied to any doubly-bounded variable to make it doubly unbounded.<br>
!>  A prime example of this is the binary mixture fraction.<br>
!>
!>  \see
!>  [pm_sampleCor](@ref pm_sampleCor)<br>
!>
!>  \test
!>  [test_pm_mathFisher](@ref test_pm_mathFisher)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_mathFisher

    use pm_kind, only: IK, RK, SK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathFisher"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the Fisher transformation of the input Fisher z value.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_mathFisher](@ref pm_mathFisher) for information about the Fisher transformation.
    !>
    !>  \param[in]  val     :   The input scalar (or array of the same shape as other array-like arguments) of,
    !>                          <ol>
    !>                              <li>    type `real` of kind \RKALL,
    !>                          </ol>
    !>                          containing the bounded value to be Fisher-transformed.
    !>  \param[in]  lb      :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `val`,
    !>                          containing the lower bound of the range of the input `val`.
    !>                          (**optional**, default = `-1`. It must be present **if and only if** the optional input argument `ub` is also present.)
    !>  \param[in]  ub      :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `val`,
    !>                          containing the upper bound of the range of the input `val`.
    !>                          (**optional**, default = `+1`. It must be present **if and only if** the optional input argument `lb` is also present.)
    !>
    !>  \return
    !>  `fisherz`           :   The output scalar (or array of the same shape as other array-like arguments) of the same type and kind as `val`,
    !>                          containing the result of the Fisher transformation of the input `val`.
    !>
    !>  \interface{getFisher}
    !>  \code{.F90}
    !>
    !>      use pm_mathFisher, only: getFisher
    !>
    !>      fisherz = getFisher(val)
    !>      fisherz = getFisher(val, lb, ub)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `lb < ub` must hold for the corresponding input arguments.<br>
    !>  The condition `lb < val .and. val < ub` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getFisher](@ref pm_mathFisher::getFisher)<br>
    !>  [getFisherInv](@ref pm_mathFisher::getFisherInv)<br>
    !>
    !>  \example{getFisher}
    !>  \include{lineno} example/pm_mathFisher/getFisher/main.F90
    !>  \compilef{getFisher}
    !>  \output{getFisher}
    !>  \include{lineno} example/pm_mathFisher/getFisher/main.out.F90
    !>  \postproc{getFisher}
    !>  \include{lineno} example/pm_mathFisher/getFisher/main.py
    !>  \vis{getFisher}
    !>  \image html pm_mathFisher/getFisher/getFisher.png width=700
    !>
    !>  \test
    !>  [test_pm_mathFisher](@ref test_pm_mathFisher)
    !>
    !>  \final{getFisher}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface getFisher

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getFisherFDD_RK5(val) result(fisherz)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFisherFDD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    :: val
        real(RKG)                   :: fisherz
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getFisherFDD_RK4(val) result(fisherz)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFisherFDD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    :: val
        real(RKG)                   :: fisherz
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getFisherFDD_RK3(val) result(fisherz)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFisherFDD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    :: val
        real(RKG)                   :: fisherz
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getFisherFDD_RK2(val) result(fisherz)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFisherFDD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    :: val
        real(RKG)                   :: fisherz
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getFisherFDD_RK1(val) result(fisherz)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFisherFDD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    :: val
        real(RKG)                   :: fisherz
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getFisherFLU_RK5(val, lb, ub) result(fisherz)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFisherFLU_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    :: val, lb, ub
        real(RKG)                   :: fisherz
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getFisherFLU_RK4(val, lb, ub) result(fisherz)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFisherFLU_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    :: val, lb, ub
        real(RKG)                   :: fisherz
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getFisherFLU_RK3(val, lb, ub) result(fisherz)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFisherFLU_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    :: val, lb, ub
        real(RKG)                   :: fisherz
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getFisherFLU_RK2(val, lb, ub) result(fisherz)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFisherFLU_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    :: val, lb, ub
        real(RKG)                   :: fisherz
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getFisherFLU_RK1(val, lb, ub) result(fisherz)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFisherFLU_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    :: val, lb, ub
        real(RKG)                   :: fisherz
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the inverse Fisher transformation of the input Fisher z value.<br>
    !>
    !>  \details
    !>  See the documentation of [pm_mathFisher](@ref pm_mathFisher) for information about the Fisher transformation.
    !>
    !>  \param[in]  fisherz :   The input scalar (or array of the same shape as other array-like arguments) of,
    !>                          <ol>
    !>                              <li>    type `real` of kind \RKALL,
    !>                          </ol>
    !>                          containing the Fisher z value(s).
    !>  \param[in]  lb      :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `fisherz`,
    !>                          containing the lower bound of the range of the output `val`.
    !>                          (**optional**, default = `-1`. It must be present **if and only if** the optional input argument `ub` is also present.)
    !>  \param[in]  ub      :   The input scalar (or array of the same shape as other array-like arguments) of the same type and kind as `fisherz`,
    !>                          containing the upper bound of the range of the output `val`.
    !>                          (**optional**, default = `+1`. It must be present **if and only if** the optional input argument `lb` is also present.)
    !>
    !>  \return
    !>  `val`               :   The output scalar (or array of the same shape as other array-like arguments) of the same type and kind as `fisherz`,
    !>                          containing the result of the inverse Fisher transformation of the input `fisherz`.
    !>
    !>  \interface{getFisherInv}
    !>  \code{.F90}
    !>
    !>      use pm_mathFisher, only: getFisherInv
    !>
    !>      val = getFisherInv(fisherz)
    !>      val = getFisherInv(fisherz, lb, ub)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `lb < ub` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>  The procedures of this generic interface are **always pure** when both optional arguments `lb` and `ub` are missing.<br>
    !>
    !>  \see
    !>  [getFisher](@ref pm_mathFisher::getFisher)<br>
    !>  [getFisherInv](@ref pm_mathFisher::getFisherInv)<br>
    !>
    !>  \example{getFisherInv}
    !>  \include{lineno} example/pm_mathFisher/getFisherInv/main.F90
    !>  \compilef{getFisherInv}
    !>  \output{getFisherInv}
    !>  \include{lineno} example/pm_mathFisher/getFisherInv/main.out.F90
    !>  \postproc{getFisherInv}
    !>  \include{lineno} example/pm_mathFisher/getFisherInv/main.py
    !>  \vis{getFisherInv}
    !>  \image html pm_mathFisher/getFisherInv/getFisherInv.png width=700
    !>
    !>  \test
    !>  [test_pm_mathFisher](@ref test_pm_mathFisher)
    !>
    !>  \final{getFisherInv}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface getFisherInv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function getFisherInvFDD_RK5(fisherz) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFisherInvFDD_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    :: fisherz
        real(RKG)                   :: val
    end function
#endif

#if RK4_ENABLED
    pure elemental module function getFisherInvFDD_RK4(fisherz) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFisherInvFDD_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    :: fisherz
        real(RKG)                   :: val
    end function
#endif

#if RK3_ENABLED
    pure elemental module function getFisherInvFDD_RK3(fisherz) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFisherInvFDD_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    :: fisherz
        real(RKG)                   :: val
    end function
#endif

#if RK2_ENABLED
    pure elemental module function getFisherInvFDD_RK2(fisherz) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFisherInvFDD_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    :: fisherz
        real(RKG)                   :: val
    end function
#endif

#if RK1_ENABLED
    pure elemental module function getFisherInvFDD_RK1(fisherz) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFisherInvFDD_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    :: fisherz
        real(RKG)                   :: val
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getFisherInvFLU_RK5(fisherz, lb, ub) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFisherInvFLU_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    :: fisherz, lb, ub
        real(RKG)                   :: val
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getFisherInvFLU_RK4(fisherz, lb, ub) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFisherInvFLU_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    :: fisherz, lb, ub
        real(RKG)                   :: val
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getFisherInvFLU_RK3(fisherz, lb, ub) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFisherInvFLU_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    :: fisherz, lb, ub
        real(RKG)                   :: val
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getFisherInvFLU_RK2(fisherz, lb, ub) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFisherInvFLU_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    :: fisherz, lb, ub
        real(RKG)                   :: val
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getFisherInvFLU_RK1(fisherz, lb, ub) result(val)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFisherInvFLU_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    :: fisherz, lb, ub
        real(RKG)                   :: val
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_mathFisher ! LCOV_EXCL_LINE