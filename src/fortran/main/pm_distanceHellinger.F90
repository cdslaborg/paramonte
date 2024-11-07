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
!>  This module contains classes and procedures for computing the Hellinger statistical distance between two probability distributions.
!>
!>  \details
!>  The Hellinger distance (which is also closely related to the Bhattacharyya distance)
!>  is used to quantify the similarity between two probability distributions.<br>
!>  The Hellinger distance is defined in terms of the Hellinger integral, which was introduced by Ernst Hellinger in 1909.<br>
!>  It is also sometimes called the Jeffreys distance.<br>
!>
!>  Definition using Measure Theory
!>  -------------------------------
!>
!>  Let \f$P\f$ and \f$Q\f$ denote two probability measures on a measure space \f$\mathcal{X}\f$
!>  that are absolutely continuous with respect to an auxiliary measure \f$\lambda\f$.<br>
!>  The square of the Hellinger distance between \f$P\f$ and \f$Q\f$ is defined as the quantity,
!>  \f{equation}{
!>      H^{2}(P, Q) = {\frac{1}{2}} \int_{\mathcal{X}} \left( {\sqrt{p(x)}} - {\sqrt {q(x)}} \right)^{2} \lambda(dx) ~.
!>  \f}
!>
!>  where, \f$P(dx) = p(x)\lambda(dx)\f$ and \f$Q(dx) = q(x)\lambda(dx)\f$, that is \f$p\f$ and \f$q\f$
!>  are the Radon–Nikodym derivatives of \f$P\f$ and \f$Q\f$ respectively with respect to \f$\lambda\f$.<br>
!>  This definition does not depend on \f$\lambda\f$, that is, the Hellinger distance between \f$P\f$ and \f$Q\f$
!>  does not change if \f$\lambda\f$ is replaced with a different probability measure with respect to which both \f$P\f$ and \f$Q\f$ are absolutely continuous.<br>
!>  For compactness, the above formula is often written as,
!>  \f{equation}{
!>      H^{2}(P,Q) = {\frac{1}{2}}\int_{\mathcal {X}}\left({\sqrt {P(dx)}}-{\sqrt {Q(dx)}}\right)^{2} ~.
!>  \f}
!>
!>  Definition using Probability Theory
!>  -----------------------------------
!>
!>  To define the Hellinger distance in terms of elementary probability theory, let \f$\lambda\f$ be the Lebesgue measure,
!>  so that \f$f = \frac{dP}{d\lambda}\f$ and \f$q = \frac{dQ}{d\lambda}\f$ are simply probability density functions.<br>
!>  The squared Hellinger distance can be then expressed as a standard calculus integral,
!>  \f{equation}{
!>      H^{2}(f,g) = {\frac {1}{2}}\int \left({\sqrt {f(x)}}-{\sqrt {g(x)}}\right)^{2}\,dx=1-\int {\sqrt {f(x)g(x)}}\,dx ~,
!>  \f}
!>
!>  where the second form can be obtained by expanding the square and using the fact that the integral of a probability density over its domain equals \f$1\f$.
!>
!>  The Hellinger distance \f$H(P, Q)\f$ satisfies the property (derivable from the Cauchy–Schwarz inequality),
!>  \f{equation}{
!>      0\leq H(P,Q)\leq 1 ~.
!>  \f}
!>
!>  Definition for Discrete distributions
!>  -------------------------------------
!>
!>  For two discrete probability distributions \f$P = (p_{1}, \ldots , p_{k})\f$ and \f$Q = (q_{1},\ldots ,q_{k})\f$, their Hellinger distance is defined as,
!>  \f{equation}{
!>      H(P,Q) = {\frac {1}{\sqrt {2}}}\;{\sqrt {\sum _{i=1}^{k}({\sqrt {p_{i}}}-{\sqrt {q_{i}}})^{2}}} ~,
!>  \f}
!>
!>  which is directly related to the Euclidean norm of the difference of the square root vectors,
!>  \f{equation}{
!>      H(P,Q) = {\frac {1}{\sqrt {2}}}\;{\bigl \|}{\sqrt {P}}-{\sqrt {Q}}{\bigr \|}_{2} ~.
!>  \f}
!>
!>  It follows that,
!>  \f{equation}{
!>      1 - H^{2}(P, Q) = \sum_{i=1}^{k}{\sqrt{p_{i} q_{i}}} ~.
!>  \f}
!>
!>  Properties of the Hellinger Distance
!>  ------------------------------------
!>
!>  <ol>
!>      <li>    The Hellinger distance forms a bounded metric on the space of probability distributions over a given probability space.<br>
!>      <li>    The maximum distance \f$1\f$ is achieved when \f$P\f$ assigns probability zero to every set to which \f$Q\f$ assigns a positive probability, and vice versa.<br>
!>      <li>    Sometimes the factor \f$\frac{1}{2}\f$ in front of the integral is omitted, in which case the Hellinger distance ranges from zero to the square root of two.<br>
!>      <li>    The Hellinger distance is related to the Bhattacharyya coefficient \f$BC(P,Q)\f$ as it can be defined as,
!>              \f{equation}{
!>                  H(P,Q) = {\sqrt{1 - BC(P,Q)}} ~.
!>              \f}
!>      <li>    Hellinger distances are used in theory of sequential and asymptotic statistics.
!>      <li>    The squared Hellinger distance between two normal distributions \f$P\sim{\mathcal{N}}(\mu_{1}, \sigma_{1}^{2})\f$ and \f$Q\sim{\mathcal{N}}(\mu_{2},\sigma_{2}^{2})\f$ is,
!>              \f{equation}{
!>                  H^{2}(P,Q) = 1 - {\sqrt {\frac {2\sigma_{1}\sigma_{2}}{\sigma_{1}^{2}+\sigma _{2}^{2}}}}\,e^{-{\frac {1}{4}}{\frac {(\mu _{1}-\mu _{2})^{2}}{\sigma _{1}^{2}+\sigma _{2}^{2}}}} ~.
!>              \f}
!>      <li>    The squared Hellinger distance between two multivariate normal distributions \f$P\sim {\mathcal {N}}(\mu _{1},\Sigma _{1})\f$ and \f$Q\sim {\mathcal {N}}(\mu _{2},\Sigma _{2})\f$ is,
!>              \f{equation}{
!>                  H^{2}(P,Q)=1-{\frac {\det(\Sigma _{1})^{1/4}\det(\Sigma _{2})^{1/4}}{\det \left({\frac {\Sigma _{1}+\Sigma _{2}}{2}}\right)^{1/2}}}\exp \left\{-{\frac {1}{8}}(\mu _{1}-\mu _{2})^{T}\left({\frac {\Sigma _{1}+\Sigma _{2}}{2}}\right)^{-1}(\mu _{1}-\mu _{2})\right\} ~.
!>              \f}
!>      <li>    The squared Hellinger distance between two exponential distributions \f$P\sim \mathrm {Exp} (\alpha )\f$ and \f$Q\sim \mathrm {Exp} (\beta )\f$ is,
!>              \f{equation}{
!>                  H^{2}(P,Q) = 1 - {\frac {2{\sqrt {\alpha \beta }}}{\alpha +\beta }} ~.
!>              \f}
!>      <li>    The squared Hellinger distance between two Weibull distributions \f$P\sim \mathrm {W} (k, \alpha)\f$ and \f$Q\sim \mathrm {W} (k, \beta)\f$,
!>              where \f$k\f$ is a common shape parameter and \f$\alpha\f$ and \f$\beta\f$ are the scale parameters respectively, is,
!>              \f{equation}{
!>                  H^{2}(P, Q) = 1 - {\frac {2(\alpha \beta )^{k/2}}{\alpha ^{k}+\beta ^{k}}} ~.
!>              \f}
!>      <li>    The squared Hellinger distance between two Poisson distributions with rate parameters \f$\alpha\f$ and \f$\beta\f$, so that \f$P\sim \mathrm{Poisson}(\alpha)\f$ and \f$Q\sim \mathrm {Poisson} (\beta)\f$, is,
!>              \f{equation}{
!>                  H^{2}(P, Q) = 1 - e^{-{\frac {1}{2}}({\sqrt {\alpha }}-{\sqrt {\beta }})^{2}} ~.
!>              \f}
!>      <li>    The squared Hellinger distance between two beta distributions \f$P\sim {\text{Beta}}(a_{1},b_{1})\f$ and \f$Q\sim {\text{Beta}}(a_{2},b_{2})\f$ is,
!>              \f{equation}{
!>                  H^{2}(P,Q) = 1 - {\frac {B\left({\frac {a_{1}+a_{2}}{2}},{\frac {b_{1}+b_{2}}{2}}\right)}{\sqrt {B(a_{1},b_{1})B(a_{2},b_{2})}}} ~,
!>              \f}
!>              where \f$B\f$ represents the [beta function](@ref pm_mathBeta).
!>      <li>    The squared Hellinger distance between two gamma distributions \f$P\sim {\text{Gamma}}(a_{1}, b_{1})\f$ and \f$Q\sim {\text{Gamma}}(a_{2},b_{2})\f$ is,
!>              \f{equation}{
!>                  H^{2}(P, Q) = 1 - \Gamma\left({\scriptstyle{\frac {a_{1}+a_{2}}{2}}}\right)\left({\frac {b_{1}+b_{2}}{2}}\right)^{-(a_{1}+a_{2})/2}{\sqrt {\frac {b_{1}^{a_{1}}b_{2}^{a_{2}}}{\Gamma (a_{1})\Gamma (a_{2})}}} ~,
!>              \f}
!>              where \f$\Gamma\f$ is the [gamma function](@ref pm_mathGamma).
!>  </ol>
!>
!>  Connection with Total Variation Distance (TVD)
!>  ----------------------------------------------
!>
!>  The Hellinger distance \f$H(P, Q)\f$ and the total variation distance (or statistical distance) \f$\delta(P,Q)\f$ are related as follows,
!>  \f{equation}{
!>      H^{2}(P, Q)\leq \delta(P, Q)\leq {\sqrt{2}}H(P, Q) ~.
!>  \f}
!>
!>  These inequalities follow immediately from the inequalities between the 1-norm and the 2-norm.
!>
!>  \see
!>  [pm_distanceBhat](@ref pm_distanceBhat)<br>
!>  [pm_distanceEuclid](@ref pm_distanceEuclid)<br>
!>  [pm_distanceHellinger](@ref pm_distanceHellinger)<br>
!>  [pm_distanceKolm](@ref pm_distanceKolm)<br>
!>  [pm_distanceMahal](@ref pm_distanceMahal)<br>
!>
!>  \test
!>  [test_pm_distanceHellinger](@ref test_pm_distanceHellinger)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, March 22, 2012, 2:21 PM, National Institute for Fusion Studies, The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distanceHellinger

    use pm_kind, only: SK, IK, LK
    use pm_mathConst, only: ninf_type
    use pm_mathConst, only: pinf_type
    use pm_except, only: getInfNeg
    use pm_except, only: getInfPos

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distanceHellinger"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the square of the Hellinger distance of two univariate (discrete or continuous) distributions.
    !>
    !>  \details
    !>  See [pm_distanceHellinger](@ref pm_distanceHellinger) for the mathematical definition of the Hellinger distance.
    !>
    !>  \param[in]  p       :   The vector of the same size as the non-zero size of the input argument `q`, of the same type and kind as the output `hellsq`,
    !>                          representing the Probability Mass Function (PMF) of the first distribution in the computation of the squared Hellinger distance.<br>
    !>                          Alternatively, the input `p` can be a non-elemental function that takes a scalar of the same type and kind as the output `hellsq`
    !>                          and returns the output `pdf` of the same type and kind as `x`, representing the Probability Density Function (PDF)
    !>                          of the first distribution in the computation of the squared Hellinger distance.<br>
    !>                          The following illustrates the generic interface of `p`,
    !>                          \code{.F90}
    !>                              function p(x) result(pdf)
    !>                                  real(RKG), intent(in) :: x
    !>                                  real(RKG) :: pdf
    !>                              end function
    !>                          \endcode
    !>                          where `RKG` is the kind type parameter of the output `hellsq`.<br>
    !>                          The input arguments `p` and `q` must be of the same type and kind (and size if they are vectors).<br>
    !>  \param[in]  q       :   The vector of the same size as the non-zero size of the input argument `p`, of the same type and kind as the output `hellsq`,
    !>                          representing the Probability Mass Function (PMF) of the second distribution in the computation of the squared Hellinger distance.<br>
    !>                          Alternatively, the input `p` can be a non-elemental function that takes a scalar of the same type and kind as the output `hellsq`
    !>                          and returns the output `pdf` of the same type and kind as `x`, representing the Probability Density Function (PDF)
    !>                          of the second distribution in the computation of the squared Hellinger distance.<br>
    !>                          The following illustrates the generic interface of `p`,
    !>                          \code{.F90}
    !>                              function p(x) result(pdf)
    !>                                  real(RKG), intent(in) :: x
    !>                                  real(RKG) :: pdf
    !>                              end function
    !>                          \endcode
    !>                          where `RKG` is the kind type parameter of the output `hellsq`.<br>
    !>                          The input arguments `p` and `q` must be of the same type and kind (and size if they are vectors).<br>
    !>  \param[in]  lb      :   The input scalar of type `real` of the same kind as `integral`, representing the lower limit of integration.<br>
    !>                          Set `lb = -huge(lb)` or to the IEEE-compliant negative infinity (`lb = `[getInfNeg(lb)](@ref pm_except::getInfNeg))
    !>                          to imply \f$-\infty\f$ as the lower bound of integration.<br>
    !>                          See also the corresponding argument of [isFailedQuad()](@ref pm_quadPack::isFailedQuad).<br>
    !>                          (**optional**, default = [getInfNeg(lb)](@ref pm_except::getInfNeg). It can be present **only if** the input arguments `p` and `q` are procedures.)
    !>  \param[in]  ub      :   The input scalar of type `real` of the same kind as `integral`, representing the upper limit of integration.<br>
    !>                          Set `ub = huge(ub)` or to the IEEE-compliant positive infinity (`ub = `[getInfPos(lb)](@ref pm_except::getInfPos))
    !>                          to imply \f$+\infty\f$ as the upper bound of integration.<br>
    !>                          See also the corresponding argument of [isFailedQuad()](@ref pm_quadPack::isFailedQuad).<br>
    !>                          (**optional**, default = [getInfPos(lb)](@ref pm_except::getInfPos). It can be present **only if** the input arguments `p` and `q` are procedures.)
    !>  \param[out] failed  :   The output scalar of type `logical` of default kind \LK, that is set to `.true.` <b>if and only if</b> the integration Hellinger **fails to converge within the tolerances**.<br>
    !>                          Otherwise, it is set `.false.` if the integration succeeds with no errors.<br>
    !>                          See also the corresponding output argument of [isFailedQuad()](@ref pm_quadPack::isFailedQuad).<br>
    !>                          See also the description of the output argument `err` of [getQuadErr](@ref pm_quadPack::getQuadErr) for information on the kinds of integration failures that can happen.<br>
    !>                          (**optional**. If missing and the integration fails, the procedure will halt the program by calling `error stop`.)<br>
    !>  \param[out] msg     :   The output scalar argument of type `character` of default kind \SK of arbitrary length type parameter that is set to a diagnostic message if the integration fails to converge.<br>
    !>                          A length type parameter of `127` is sufficient to capture all error messages.<br>
    !>                          If `msg` has shorter length parameter, the output message will be trimmed from the end, otherwise padded with blanks as necessary.<br>
    !>                          See also the corresponding argument of [isFailedQuad()](@ref pm_quadPack::isFailedQuad).<br>
    !>                          (**optional**. If missing, no diagnostic message will be returned.)<br>
    !>
    !>  \return
    !>  `hellsq`            :   The output scalar of<br>
    !>                          <ol>
    !>                              <li>    type `real` of kind \RKALL,<br>
    !>                          </ol>
    !>                          representing the square of the Hellinger distance of the discrete or continuous distributions represented by `p` and `q`.<br>
    !>
    !>  \interface{getDisHellSq}
    !>  \code{.F90}
    !>
    !>      use pm_distanceHellinger, only: getDisHellSq
    !>
    !>      ! discrete distributions.
    !>
    !>      hellsq = getDisHellSq(p(:), q(:))
    !>
    !>      ! continuous univariate distributions.
    !>
    !>      hellsq = getDisHellSq(p, q, lb = lb, ub = ub)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `sum(p) == 1` must hold for the corresponding input arguments.<br>
    !>  The condition `sum(q) == 1` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0 <= p) .and. all(p <= 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0 <= q) .and. all(q <= 1)` must hold for the corresponding input arguments.<br>
    !>  The condition `size(p) == size(q)` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>  The procedures under this generic interface are always `impure` when the input arguments `p` and `q` are procedures.
    !>
    !>  \example{getDisHellSq}
    !>  \include{lineno} example/pm_distanceHellinger/getDisHellSq/main.F90
    !>  \compilef{getDisHellSq}
    !>  \output{getDisHellSq}
    !>  \include{lineno} example/pm_distanceHellinger/getDisHellSq/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distanceHellinger](@ref test_pm_distanceHellinger)<br>
    !>
    !>  \todo
    !>  \pvhigh
    !>  The runtime checks for the complex input `invCov` must be implemented.<br>
    !>
    !>  \final{getDisHellSq}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

    ! PMF

    interface getDisHellSq

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDisHellSqPMF_RK5(p, q) result(hellsq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPMF_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in), contiguous    :: p(:), q(:)
        real(RKG)                               :: hellsq
    end function
#endif

#if RK4_ENABLED
    PURE module function getDisHellSqPMF_RK4(p, q) result(hellsq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPMF_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in), contiguous    :: p(:), q(:)
        real(RKG)                               :: hellsq
    end function
#endif

#if RK3_ENABLED
    PURE module function getDisHellSqPMF_RK3(p, q) result(hellsq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPMF_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in), contiguous    :: p(:), q(:)
        real(RKG)                               :: hellsq
    end function
#endif

#if RK2_ENABLED
    PURE module function getDisHellSqPMF_RK2(p, q) result(hellsq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPMF_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in), contiguous    :: p(:), q(:)
        real(RKG)                               :: hellsq
    end function
#endif

#if RK1_ENABLED
    PURE module function getDisHellSqPMF_RK1(p, q) result(hellsq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPMF_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in), contiguous    :: p(:), q(:)
        real(RKG)                               :: hellsq
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! Proc

    interface getDisHellSq

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getDisHellSqPRC_RK5(p, q, lb, ub, failed, msg) result(hellsq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPRC_RK5
#endif
        use pm_kind, only: RKG => RK5
        procedure(real(RKG))                        :: p, q
        real(RKG)       , intent(in)    , optional  :: lb
        real(RKG)       , intent(in)    , optional  :: ub
        character(*, SK), intent(out)   , optional  :: msg
        logical(LK)     , intent(out)   , optional  :: failed
        real(RKG)                                   :: hellsq
    end function
#endif

#if RK4_ENABLED
    module function getDisHellSqPRC_RK4(p, q, lb, ub, failed, msg) result(hellsq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPRC_RK4
#endif
        use pm_kind, only: RKG => RK4
        procedure(real(RKG))                        :: p, q
        real(RKG)       , intent(in)    , optional  :: lb
        real(RKG)       , intent(in)    , optional  :: ub
        character(*, SK), intent(out)   , optional  :: msg
        logical(LK)     , intent(out)   , optional  :: failed
        real(RKG)                                   :: hellsq
    end function
#endif

#if RK3_ENABLED
    module function getDisHellSqPRC_RK3(p, q, lb, ub, failed, msg) result(hellsq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPRC_RK3
#endif
        use pm_kind, only: RKG => RK3
        procedure(real(RKG))                        :: p, q
        real(RKG)       , intent(in)    , optional  :: lb
        real(RKG)       , intent(in)    , optional  :: ub
        character(*, SK), intent(out)   , optional  :: msg
        logical(LK)     , intent(out)   , optional  :: failed
        real(RKG)                                   :: hellsq
    end function
#endif

#if RK2_ENABLED
    module function getDisHellSqPRC_RK2(p, q, lb, ub, failed, msg) result(hellsq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPRC_RK2
#endif
        use pm_kind, only: RKG => RK2
        procedure(real(RKG))                        :: p, q
        real(RKG)       , intent(in)    , optional  :: lb
        real(RKG)       , intent(in)    , optional  :: ub
        character(*, SK), intent(out)   , optional  :: msg
        logical(LK)     , intent(out)   , optional  :: failed
        real(RKG)                                   :: hellsq
    end function
#endif

#if RK1_ENABLED
    module function getDisHellSqPRC_RK1(p, q, lb, ub, failed, msg) result(hellsq)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPRC_RK1
#endif
        use pm_kind, only: RKG => RK1
        procedure(real(RKG))                        :: p, q
        real(RKG)       , intent(in)    , optional  :: lb
        real(RKG)       , intent(in)    , optional  :: ub
        character(*, SK), intent(out)   , optional  :: msg
        logical(LK)     , intent(out)   , optional  :: failed
        real(RKG)                                   :: hellsq
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    module function getDisHellSqPRC_FF_RK5(p, q, lb, ub) result(hellsq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPRC_FF_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        procedure(real(RKG))                    :: p, q
!        real(RKG)       , intent(in)            :: lb
!        real(RKG)       , intent(in)            :: ub
!        real(RKG)                               :: hellsq
!    end function
!#endif
!
!#if RK4_ENABLED
!    module function getDisHellSqPRC_FF_RK4(p, q, lb, ub) result(hellsq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPRC_FF_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        procedure(real(RKG))                    :: p, q
!        real(RKG)       , intent(in)            :: lb
!        real(RKG)       , intent(in)            :: ub
!        real(RKG)                               :: hellsq
!    end function
!#endif
!
!#if RK3_ENABLED
!    module function getDisHellSqPRC_FF_RK3(p, q, lb, ub) result(hellsq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPRC_FF_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        procedure(real(RKG))                    :: p, q
!        real(RKG)       , intent(in)            :: lb
!        real(RKG)       , intent(in)            :: ub
!        real(RKG)                               :: hellsq
!    end function
!#endif
!
!#if RK2_ENABLED
!    module function getDisHellSqPRC_FF_RK2(p, q, lb, ub) result(hellsq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPRC_FF_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        procedure(real(RKG))                    :: p, q
!        real(RKG)       , intent(in)            :: lb
!        real(RKG)       , intent(in)            :: ub
!        real(RKG)                               :: hellsq
!    end function
!#endif
!
!#if RK1_ENABLED
!    module function getDisHellSqPRC_FF_RK1(p, q, lb, ub) result(hellsq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPRC_FF_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        procedure(real(RKG))                    :: p, q
!        real(RKG)       , intent(in)            :: lb
!        real(RKG)       , intent(in)            :: ub
!        real(RKG)                               :: hellsq
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    module function getDisHellSqPRC_FI_RK5(p, q, lb, ub) result(hellsq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPRC_FI_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        procedure(real(RKG))                    :: p, q
!        real(RKG)       , intent(in)            :: lb
!        type(pinf_type) , intent(in)            :: ub
!        real(RKG)                               :: hellsq
!    end function
!#endif
!
!#if RK4_ENABLED
!    module function getDisHellSqPRC_FI_RK4(p, q, lb, ub) result(hellsq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPRC_FI_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        procedure(real(RKG))                    :: p, q
!        real(RKG)       , intent(in)            :: lb
!        type(pinf_type) , intent(in)            :: ub
!        real(RKG)                               :: hellsq
!    end function
!#endif
!
!#if RK3_ENABLED
!    module function getDisHellSqPRC_FI_RK3(p, q, lb, ub) result(hellsq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPRC_FI_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        procedure(real(RKG))                    :: p, q
!        real(RKG)       , intent(in)            :: lb
!        type(pinf_type) , intent(in)            :: ub
!        real(RKG)                               :: hellsq
!    end function
!#endif
!
!#if RK2_ENABLED
!    module function getDisHellSqPRC_FI_RK2(p, q, lb, ub) result(hellsq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPRC_FI_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        procedure(real(RKG))                    :: p, q
!        real(RKG)       , intent(in)            :: lb
!        type(pinf_type) , intent(in)            :: ub
!        real(RKG)                               :: hellsq
!    end function
!#endif
!
!#if RK1_ENABLED
!    module function getDisHellSqPRC_FI_RK1(p, q, lb, ub) result(hellsq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPRC_FI_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        procedure(real(RKG))                    :: p, q
!        real(RKG)       , intent(in)            :: lb
!        type(pinf_type) , intent(in)            :: ub
!        real(RKG)                               :: hellsq
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    module function getDisHellSqPRC_IF_RK5(p, q, lb, ub) result(hellsq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPRC_IF_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        procedure(real(RKG))                    :: p, q
!        type(ninf_type) , intent(in)            :: lb
!        real(RKG)       , intent(in)            :: ub
!        real(RKG)                               :: hellsq
!    end function
!#endif
!
!#if RK4_ENABLED
!    module function getDisHellSqPRC_IF_RK4(p, q, lb, ub) result(hellsq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPRC_IF_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        procedure(real(RKG))                    :: p, q
!        type(ninf_type) , intent(in)            :: lb
!        real(RKG)       , intent(in)            :: ub
!        real(RKG)                               :: hellsq
!    end function
!#endif
!
!#if RK3_ENABLED
!    module function getDisHellSqPRC_IF_RK3(p, q, lb, ub) result(hellsq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPRC_IF_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        procedure(real(RKG))                    :: p, q
!        type(ninf_type) , intent(in)            :: lb
!        real(RKG)       , intent(in)            :: ub
!        real(RKG)                               :: hellsq
!    end function
!#endif
!
!#if RK2_ENABLED
!    module function getDisHellSqPRC_IF_RK2(p, q, lb, ub) result(hellsq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPRC_IF_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        procedure(real(RKG))                    :: p, q
!        type(ninf_type) , intent(in)            :: lb
!        real(RKG)       , intent(in)            :: ub
!        real(RKG)                               :: hellsq
!    end function
!#endif
!
!#if RK1_ENABLED
!    module function getDisHellSqPRC_IF_RK1(p, q, lb, ub) result(hellsq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPRC_IF_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        procedure(real(RKG))                    :: p, q
!        type(ninf_type) , intent(in)            :: lb
!        real(RKG)       , intent(in)            :: ub
!        real(RKG)                               :: hellsq
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    module function getDisHellSqPRC_II_RK5(p, q, lb, ub) result(hellsq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPRC_II_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        procedure(real(RKG))                    :: p, q
!        type(ninf_type) , intent(in)            :: ninf_type
!        type(pinf_type) , intent(in)            :: pinf_type
!        real(RKG)                               :: hellsq
!    end function
!#endif
!
!#if RK4_ENABLED
!    module function getDisHellSqPRC_II_RK4(p, q, lb, ub) result(hellsq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPRC_II_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        procedure(real(RKG))                    :: p, q
!        type(ninf_type) , intent(in)            :: ninf_type
!        type(pinf_type) , intent(in)            :: pinf_type
!        real(RKG)                               :: hellsq
!    end function
!#endif
!
!#if RK3_ENABLED
!    module function getDisHellSqPRC_II_RK3(p, q, lb, ub) result(hellsq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPRC_II_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        procedure(real(RKG))                    :: p, q
!        type(ninf_type) , intent(in)            :: ninf_type
!        type(pinf_type) , intent(in)            :: pinf_type
!        real(RKG)                               :: hellsq
!    end function
!#endif
!
!#if RK2_ENABLED
!    module function getDisHellSqPRC_II_RK2(p, q, lb, ub) result(hellsq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPRC_II_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        procedure(real(RKG))                    :: p, q
!        type(ninf_type) , intent(in)            :: ninf_type
!        type(pinf_type) , intent(in)            :: pinf_type
!        real(RKG)                               :: hellsq
!    end function
!#endif
!
!#if RK1_ENABLED
!    module function getDisHellSqPRC_II_RK1(p, q, lb, ub) result(hellsq)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getDisHellSqPRC_II_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        procedure(real(RKG))                    :: p, q
!        type(ninf_type) , intent(in)            :: ninf_type
!        type(pinf_type) , intent(in)            :: pinf_type
!        real(RKG)                               :: hellsq
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distanceHellinger ! LCOV_EXCL_LINE